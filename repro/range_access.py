"""Budget-limited, auditable HTTP Range access to public Figshare HDF5 files.

Reads only requested byte blocks; refuses a server returning an entire file.
"""
import io
import json
import pathlib
import time
import requests
import h5py


class RangeFile(io.RawIOBase):
    def __init__(self, file_id, size, cache_dir, budget=490_000_000, block_size=65536):
        self.file_id, self.size, self.pos = int(file_id), int(size), 0
        self.block_size, self.budget, self.transferred = block_size, budget, 0
        self.url = f"https://ndownloader.figshare.com/files/{file_id}"
        self.cache = pathlib.Path(cache_dir) / str(file_id)
        self.cache.mkdir(parents=True, exist_ok=True)
        self.session = requests.Session()
        self.requests = 0

    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.pos
    def seek(self, offset, whence=0):
        self.pos = offset if whence == 0 else self.pos + offset if whence == 1 else self.size + offset
        return self.pos

    def _block(self, idx):
        path = self.cache / f"{idx}.bin"
        if path.exists(): return path.read_bytes()
        start, end = idx * self.block_size, min(self.size, (idx + 1) * self.block_size) - 1
        length = end - start + 1
        if self.transferred + length > self.budget:
            raise RuntimeError("Transfer budget reached")
        for attempt in range(4):
            response = self.session.get(self.url, headers={"Range": f"bytes={start}-{end}"}, timeout=60, stream=True)
            if response.status_code == 206:
                break
            response.close()
            if attempt == 3: raise RuntimeError(f"Expected HTTP 206, got {response.status_code}")
            time.sleep(attempt + 1)
        content_range = response.headers.get("Content-Range", "")
        if not content_range.startswith(f"bytes {start}-{end}/"):
            response.close()
            raise RuntimeError(f"Incorrect range: {content_range}")
        data = response.content
        if len(data) != length: raise RuntimeError("Incomplete byte range")
        path.write_bytes(data)
        self.transferred += length
        self.requests += 1
        return data

    def read(self, n=-1):
        if n < 0: n = self.size - self.pos
        stop = min(self.size, self.pos + n)
        result = bytearray()
        while self.pos < stop:
            idx, offset = divmod(self.pos, self.block_size)
            block = self._block(idx)
            take = min(stop - self.pos, len(block) - offset)
            result.extend(block[offset:offset + take])
            self.pos += take
        return bytes(result)

    def readinto(self, b):
        data = self.read(len(b))
        b[:len(data)] = data
        return len(data)


def inspect_file(file_id, size, cache_dir):
    remote = RangeFile(file_id, size, cache_dir, budget=25_000_000)
    result = {}
    with h5py.File(remote, "r") as handle:
        def collect(name, obj):
            if isinstance(obj, h5py.Dataset):
                result[name] = {"shape": obj.shape, "dtype": str(obj.dtype), "chunks": obj.chunks,
                                "compression": obj.compression, "offset": obj.id.get_offset()}
                if obj.size <= 10:
                    result[name]["value"] = obj[()].tolist()
        handle.visititems(collect)
    result["_network"] = {"bytes": remote.transferred, "requests": remote.requests}
    return result


if __name__ == "__main__":
    import sys
    result = inspect_file(int(sys.argv[1]), int(sys.argv[2]), pathlib.Path(__file__).parent / "range-cache")
    print(json.dumps(result, indent=2), flush=True)
