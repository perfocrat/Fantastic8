# Fantastic8: Quantum Well Explorer

Fantastic8 calculates states for two-dimensional infinite-potential wells in six geometries. The HTML, CSS, and JavaScript interface is served by a small Python standard-library HTTP server and calls the existing shape adapters through a same-origin JSON API.

## Requirements

- Python 3.10 or newer
- Packages listed in `requirements.txt` (NumPy, SciPy, Matplotlib, and SymPy)
- Node.js for the optional JavaScript syntax check

Install and start from the project root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt
python3 server.py
```

Open <http://127.0.0.1:8000>. The server binds to localhost by default. Use `python3 server.py --host 127.0.0.1 --port 8001` to choose another local port.

## Architecture

`server.py` serves `frontend/index.html`, `style.css`, and `script.js`, and exposes `GET /api/health` and `POST /api/calculate`. The browser submits one geometry, its dimensions in nanometres, and the selected quantum numbers. The API validates the request, converts lengths to metres, and calls `shapes.adapters.calculate_shape` once; that dispatcher invokes only the selected shape adapter. The adapters and mathematical shape modules remain the physics source of truth.

The response includes the selected energy in joules and meV, available state energies, the wavefunction, probability density, coordinates, and a potential-domain grid. Potential is zero inside the well; `null` grid entries denote the infinite wall outside. The browser renders the returned grids directly and does not substitute demo results.

Geometry fields are shape-specific: rectangle uses width and height; triangle uses side; circle uses radius; sector uses radius, angle, and starting angle; annulus uses inner and outer radius; ellipse uses semi-major and semi-minor axes. Lengths are entered in nm. Quantum-number conventions follow the adapters: triangle requires `1 <= n <= m <= 5`; circle, annulus, and ellipse allow `m >= 0, n >= 1`; rectangle and sector use positive indices. Annulus maps these fields to angular momentum and radial index.

The HTTP server is deliberately single-threaded because the existing ellipse and sector adapters temporarily set module-level solver parameters. This serializes calculations so those settings cannot overlap between requests.

## Tests

Run the Python integration tests and optional JavaScript syntax check from the project root:

```bash
python3 -m unittest discover -s tests -v
node --check frontend/script.js
```

Generated physics plots are written under `generated_plots/` and are ignored by Git.