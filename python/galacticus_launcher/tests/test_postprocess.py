"""Tests for the launcher's `postprocess` module and subcommand.

These build a small stand-in for a Galacticus output file (a `Parameters` group, a `Version` group, and one output of
`nodeData`), and check the parameter file constructed for postprocessing and the validation of galaxy selections. Running
Galacticus itself is covered by `testSuite/test-postprocess-outputs.py`.
"""

import h5py
import lxml.etree as ET
import numpy as np
import pytest

from galacticus_launcher import cli
from galacticus_launcher import postprocess
from Galacticus.Parameters.extract import parameters_from_hdf5

EXTRACTOR = '<nodePropertyExtractor value="massStellar"/>'


@pytest.fixture
def catalog(tmp_path):
    """A minimal Galacticus output file with two outputs, of 5 and 3 galaxies."""
    path = tmp_path / "catalog.hdf5"
    with h5py.File(path, "w") as file:
        file.create_group("Version").attrs["gitHash"] = np.bytes_(b"0123456789abcdef0123456789abcdef01234567")
        parameters = file.create_group("Parameters")
        parameters.attrs["task"] = np.bytes_(b"evolveForests")
        parameters.attrs["galacticFilter"] = np.bytes_(b"lightcone")
        parameters.attrs["outputFileName"] = np.bytes_(b"original.hdf5")
        parameters.attrs["mergerTreeOutputter"] = np.bytes_(b"null")
        parameters.attrs["nodePropertyExtractor"] = np.bytes_(b"nodeIndices")
        parameters.attrs["starFormationHistory"] = np.bytes_(b"fixedAges")
        history = parameters.create_group("starFormationHistory")
        history.attrs["countAges"] = np.int32(10)
        # A parameter which took its default value.
        history.attrs["ageMinimum"] = np.float64(0.001)
        history.attrs["ageMinimum{defaulted}"] = np.int32(1)
        # An object which took its default value, with a sub-parameter.
        parameters.attrs["accretionHaloTotal"] = np.bytes_(b"simple")
        parameters.attrs["accretionHaloTotal{defaulted}"] = np.int32(1)
        parameters.create_group("accretionHaloTotal").attrs["timeReionization"] = np.float64(0.5)
        for index, count in ((1, 5), (2, 3)):
            file.create_dataset(f"Lightcone/Output{index}/nodeData/basicMass", data=np.ones(count))
    return path


def test_build_parameters(catalog, tmp_path):
    """The task, outputter, extractors, filter, and output file are replaced; the rest of the model is kept."""
    parameters = postprocess.build_parameters(catalog, EXTRACTOR, tmp_path / "out.hdf5", outputs_group="Lightcone",
                                              selection_file=tmp_path / "selection.hdf5", count_batch=64)
    assert len(parameters.findall("task")) == 1
    task = parameters.find("task")
    assert task.get("value") == "postprocessOutputs"
    assert task.find("fileName").get("value") == str(catalog.resolve())
    assert task.find("outputsGroupName").get("value") == "Lightcone"
    assert task.find("selectionFileName").get("value") == str((tmp_path / "selection.hdf5").resolve())
    assert task.find("countBatch").get("value") == "64"
    outputter = parameters.find("mergerTreeOutputter")
    assert outputter.get("value") == "standard"
    assert outputter.find("outputComponentProperties").get("value") == "false"
    assert outputter.find("outputsGroupName").get("value") == "Lightcone"
    extractors = parameters.findall("nodePropertyExtractor")
    assert len(extractors) == 1
    assert extractors[0].get("value") == "multi"
    assert [child.get("value") for child in extractors[0]] == ["indexPostprocessSource", "massStellar"]
    assert parameters.find("galacticFilter").get("value") == "always"
    assert parameters.find("outputFileName").get("value") == str((tmp_path / "out.hdf5").resolve())
    # The model's star formation history is kept, with defaulted parameters omitted.
    history = parameters.find("starFormationHistory")
    assert history.get("value") == "fixedAges"
    assert history.find("countAges").get("value") == "10"
    assert history.find("ageMinimum") is None
    assert parameters.find("accretionHaloTotal") is None
    assert parameters.find("lastModified").get("revision") == "0123456789abcdef0123456789abcdef01234567"


def test_build_parameters_options(catalog, tmp_path):
    """No selection, several extractors (given as elements or a file), and the filter kept or replaced."""
    extractorFile = tmp_path / "extractor.xml"
    extractorFile.write_text('<nodePropertyExtractor value="inclination"/>')
    parameters = postprocess.build_parameters(catalog, [ET.fromstring(EXTRACTOR), extractorFile], tmp_path / "out.hdf5",
                                              galactic_filter="keep")
    assert parameters.find("task/selectionFileName") is None
    assert [child.get("value") for child in parameters.find("nodePropertyExtractor")] == \
        ["indexPostprocessSource", "massStellar", "inclination"]
    assert parameters.find("galacticFilter").get("value") == "lightcone"
    parameters = postprocess.build_parameters(catalog, EXTRACTOR, tmp_path / "out.hdf5",
                                              galactic_filter='<galacticFilter value="not"/>')
    assert parameters.find("galacticFilter").get("value") == "not"


def test_build_parameters_rejects_wrong_element(catalog, tmp_path):
    with pytest.raises(ValueError, match="nodePropertyExtractor"):
        postprocess.build_parameters(catalog, '<galacticFilter value="always"/>', tmp_path / "out.hdf5")
    with pytest.raises(ValueError, match="at least one"):
        postprocess.build_parameters(catalog, [], tmp_path / "out.hdf5")


def test_extract_includes_defaulted_by_default(catalog):
    """The extraction used by `parametersExtract.py` still includes defaulted parameters unless asked not to."""
    parameters = parameters_from_hdf5(str(catalog))
    assert parameters.find("starFormationHistory/ageMinimum").get("value") == "0.001"
    assert parameters.find("accretionHaloTotal/timeReionization").get("value") == "0.5"


def test_write_selection(catalog, tmp_path):
    path = postprocess.write_selection(tmp_path / "selection.hdf5", {1: [4, 0, 4], 2: np.array([2])},
                                       outputs_group="Lightcone", catalog=catalog)
    with h5py.File(path, "r") as file:
        assert list(file["Lightcone/Output1/selection"][:]) == [4, 0, 4]
        assert list(file["Lightcone/Output2/selection"][:]) == [2]
        assert file["Lightcone/Output1/selection"].dtype == np.int64


@pytest.mark.parametrize("selection, message", [
    ({1: [5]}, "range"),
    ({1: [-1]}, "range"),
    ({3: [0]}, "no output"),
    ({1: [[0, 1]]}, "one-dimensional"),
    ({1: [0.5]}, "integers"),
    ([0, 1], "dictionary"),
])
def test_write_selection_invalid(catalog, tmp_path, selection, message):
    with pytest.raises(ValueError, match=message):
        postprocess.write_selection(tmp_path / "selection.hdf5", selection, outputs_group="Lightcone", catalog=catalog)


def test_write_selection_array_single_output(tmp_path):
    """An array alone is accepted where the file has a single output."""
    path = tmp_path / "single.hdf5"
    with h5py.File(path, "w") as file:
        file.create_dataset("Outputs/Output7/nodeData/basicMass", data=np.ones(4))
    postprocess.write_selection(tmp_path / "selection.hdf5", np.array([3, 1]), catalog=path)
    with h5py.File(tmp_path / "selection.hdf5", "r") as file:
        assert list(file["Outputs/Output7/selection"][:]) == [3, 1]


def test_cli_postprocess_dispatch(monkeypatch, tmp_path):
    """`galacticus postprocess` passes its options through to `postprocess.postprocess`."""
    calls = {}
    monkeypatch.setattr(cli, "_ensure", lambda install, tools=None: None)
    monkeypatch.setattr(cli.paths, "resolve", lambda: "install")

    def fake(catalog, extractors, output, **options):
        calls.update(catalog=catalog, extractors=extractors, output=output, **options)
        return output
    monkeypatch.setattr(postprocess, "postprocess", fake)
    status = cli.main(["postprocess", "catalog.hdf5", "extractors.xml", "change.xml", "-o", "out.hdf5",
                       "--outputs-group", "Lightcone", "--threads", "4", "--keep-filter", "--selection", "sel.hdf5"])
    assert status == 0
    assert calls["catalog"] == "catalog.hdf5"
    assert calls["extractors"] == "extractors.xml"
    assert calls["change_files"] == ["change.xml"]
    assert calls["outputs_group"] == "Lightcone"
    assert calls["threads"] == 4
    assert calls["galactic_filter"] == "keep"
    assert calls["selection"] == "sel.hdf5"
    assert calls["install"] == "install"


def test_build_parameters_dereferences(catalog, tmp_path):
    """A reference whose target lies within a replaced parameter is replaced by a copy of that target, not left dangling."""
    with h5py.File(catalog, "a") as file:
        parameters = file["Parameters"]
        parameters.create_group("task").attrs["cosmologyParameters{id:cosmology}"] = np.bytes_(b"simple")
        parameters["task"].create_group("cosmologyParameters").attrs["HubbleConstant"] = np.float64(70.0)
        parameters.attrs["cosmologyParameters"] = np.bytes_(b"{idRef:cosmology}")
    parameters = postprocess.build_parameters(catalog, EXTRACTOR, tmp_path / "out.hdf5")
    assert parameters.find("task").get("value") == "postprocessOutputs"
    cosmology = parameters.find("cosmologyParameters")
    assert cosmology.get("idRef") is None
    assert cosmology.get("value") == "simple"
    assert cosmology.find("HubbleConstant").get("value") == "70.0"
    assert parameters.findall(".//*[@idRef]") == []


def test_build_parameters_extractor_references(catalog, tmp_path):
    """Extractors may refer to parameters of the model, and to each other; the caller's elements are not altered."""
    with h5py.File(catalog, "a") as file:
        parameters = file["Parameters"]
        parameters.attrs["cosmologyFunctions{id:functions}"] = np.bytes_(b"matterLambda")
    extractors = ET.fromstring(
        '<list>'
        '<nodePropertyExtractor value="lightcone"><cosmologyFunctions idRef="functions"/></nodePropertyExtractor>'
        '<nodePropertyExtractor value="dustAttenuation">'
        '<dustAttenuation id="dust" value="atlasCompendium"><fileName value="table.hdf5"/></dustAttenuation>'
        '</nodePropertyExtractor>'
        '<nodePropertyExtractor value="dustAttenuation"><dustAttenuation idRef="dust"/></nodePropertyExtractor>'
        '</list>'
    )
    given = list(extractors)
    parameters = postprocess.build_parameters(catalog, given, tmp_path / "out.hdf5")
    assert parameters.findall(".//*[@idRef]") == []
    assert parameters.findall(".//*[@id]") == []
    multi = parameters.find("nodePropertyExtractor")
    assert multi[1].find("cosmologyFunctions").get("value") == "matterLambda"
    assert multi[3].find("dustAttenuation").get("value") == "atlasCompendium"
    assert multi[3].find("dustAttenuation/fileName").get("value") == "table.hdf5"
    # The model's own parameter is kept, without its `id`.
    assert parameters.find("cosmologyFunctions").get("value") == "matterLambda"
    # The caller's elements are unchanged, and still in their own tree.
    assert [element.getparent() for element in given] == [extractors] * 3
    assert given[0].find("cosmologyFunctions").get("idRef") == "functions"


def test_dereference_parameters_nested_and_cycles():
    """References inside a copied parameter are resolved; a cycle of references is reported."""
    from Galacticus.Parameters.extract import dereference_parameters
    parameters = ET.fromstring(
        '<parameters>'
        '<a id="x" value="1"><b idRef="y"/></a>'
        '<b id="y" value="2"/>'
        '<c><a idRef="x"/></c>'
        '</parameters>'
    )
    dereference_parameters(parameters)
    assert parameters.find("c/a").get("value") == "1"
    assert parameters.find("c/a/b").get("value") == "2"
    assert parameters.findall(".//*[@idRef]") == [] and parameters.findall(".//*[@id]") == []
    cycle = ET.fromstring('<parameters><a id="x" value="1"><a idRef="x"/></a></parameters>')
    with pytest.raises(ValueError, match="cycle"):
        dereference_parameters(cycle)
