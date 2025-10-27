# tests/conftest.py
import pytest
import shutil
import os
from pathlib import Path

@pytest.fixture
def test_env(tmp_path):
    """
    Sets up a per-test working directory with required data and config paths.
    Returns a dictionary with paths to use in tests.
    """
    root = Path(__file__).resolve().parents[1]  # project root (PhylUp)

    # Create dedicated working folder
    workdir = tmp_path / "workdir"
    workdir.mkdir()

    # tmp folder inside workdir
    tmp_folder = workdir / "tmp"
    tmp_folder.mkdir()

    # Copy template data
    shutil.copytree(root / "data/tmp_for_test", tmp_folder, dirs_exist_ok=True)

    # Paths to config and test datasets
    paths = {
        "workdir": workdir,
        "tmp_folder": tmp_folder,
        "configfi": root / "data/localblast_test.config",
        "trfn": root / "data/tiny_test_example/test.tre",
        "schema_trf": "newick",
        "id_to_spn": root / "data/tiny_test_example/test_nicespl.csv",
        "seqaln": root / "data/tiny_test_example/test.fas",
        "blast_folder": root / "data/blast_for_tests"
    }
    return paths



@pytest.fixture(autouse=True)
def set_cwd_to_project_root(monkeypatch):
    project_root = Path(__file__).parent.parent
    monkeypatch.chdir(project_root)

def pytest_configure():
    os.environ["PATH"] = "/home/blubb/PhylogeneticSoftware/PhylUp/.venv3/bin:/home/blubb/PhylogeneticSoftware/PhylUp/.venv3/bin:/home/linuxbrew/.linuxbrew/bin:/home/linuxbrew/.linuxbrew/sbin:/home/blubb/.local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/snap/bin:/home/blubb/PhylogeneticSoftware/modeltest-ng/modeltest/build:/home/blubb/PhylogeneticSoftware/PhylUp/PaPaRa:/home/blubb/PhylogeneticSoftware/EPA-ng/epa-ng-master/bin:/home/blubb/PhylogeneticSoftware/RAxML-ng/bin:/home/blubb/PhylogeneticSoftware/PhylUp/PaPaRa"
