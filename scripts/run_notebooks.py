"""Execute and save all project notebooks in order from the repository root."""
from pathlib import Path
import nbformat
from nbclient import NotebookClient

root=Path.cwd()
for path in sorted(Path("notebooks").glob("0*.ipynb")):
    print("Executing",path,flush=True)
    notebook=nbformat.read(path,as_version=4)
    NotebookClient(notebook,timeout=180,kernel_name="python3",
                   resources={"metadata":{"path":str(root)}}).execute()
    nbformat.write(notebook,path)
    print("Passed",path,flush=True)
