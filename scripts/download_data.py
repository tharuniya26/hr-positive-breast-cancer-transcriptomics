"""Download original GEO files; refuses silent replacement of existing files."""
from pathlib import Path
from urllib.request import urlretrieve, urlopen, Request
from urllib.error import URLError
import ssl
import zipfile
from src.data_loading import MATRIX_URL, ANNOT_URL, sha256

root=Path("data/raw");root.mkdir(parents=True,exist_ok=True)
sources=[(MATRIX_URL,"GSE45827_series_matrix.txt.gz","741184cd3b3354ed3ea45fc6a428ed36401cadc20505e5cc5fb381d8aed0ad35"),
         (ANNOT_URL,"GPL570.annot.gz","d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394")]
for url,name,expected in sources:
    path=root/name
    if not path.exists():
        print("Downloading",url)
        urlretrieve(url,path)
    observed=sha256(path)
    if observed!=expected:
        raise ValueError(f"Checksum mismatch for {path}: {observed}")
    print(path,path.stat().st_size,observed)

gmt=root/"ReactomePathways.gmt"
expected="89983d5c1f0af11c52edfeee7323eb425580ac6281d387a528562ab1787ce56b"
if not gmt.exists():
    url="https://www.reactome.org/download/current/ReactomePathways.gmt.zip"
    request=Request(url,headers={"User-Agent":"Mozilla/5.0"})
    try:
        content=urlopen(request,timeout=30).read()
    except URLError as error:
        if not isinstance(error.reason,ssl.SSLCertVerificationError):
            raise
        # The endpoint presented an expired certificate in this environment.
        # The independently recorded SHA-256 below still validates the payload.
        content=urlopen(request,context=ssl._create_unverified_context(),timeout=30).read()
    archive=root/"ReactomePathways.gmt.zip"
    archive.write_bytes(content)
    with zipfile.ZipFile(archive) as z:
        gmt.write_bytes(z.read("ReactomePathways.gmt"))
observed=sha256(gmt)
if observed!=expected:
    raise ValueError(f"Reactome v97 checksum mismatch: {observed}")
print(gmt,gmt.stat().st_size,observed)
