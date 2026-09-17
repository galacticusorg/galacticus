# Galacticus model tutorials

Jupyter notebooks which analyze the output of Galacticus models, using the
[dendros](https://github.com/galacticusorg/dendros) package. Each notebook is committed with executed output, so you can
read the results without running anything.

| Notebook | Model | What it covers |
|---|---|---|
| [01 — Dust absorption and emission](01-dust-absorption-and-emission.ipynb) | `parameters/tutorials/dustEmission.xml` | The spectrum of a galaxy's stars, emission lines, and AGN before and after dust; the light absorbed by birth clouds and the diffuse interstellar medium; its re-emission in the infrared by PAHs and larger grains; and energy balance between the two. |

## Running a tutorial

1. Build Galacticus (`make Galacticus.exe`), and clone the
   [datasets repository](https://github.com/galacticusorg/datasets), pointing `GALACTICUS_DATA_PATH` at it.
2. From the root of the Galacticus source tree, run the model named in the table, e.g.

   ```sh
   ./Galacticus.exe parameters/tutorials/dustEmission.xml
   ```

   which writes its output (here `dustEmissionTutorial.hdf5`) to that directory. The first run can take tens of
   minutes while Galacticus tabulates stellar population spectra and emission lines for the model; these are cached,
   so later runs take seconds.
3. Install the Python requirements: `pip install 'dendros[plot]' jupyter`.
4. Start Jupyter from this `tutorials/models/` directory and open the notebook. It looks for the model output in the
   root of the source tree; set `GALACTICUS_EXEC_PATH` to that directory if you run it from elsewhere.
