# Electric Field Visualization

This is an interactive 3D simulation of the electric field from a charged ring or disk, emphasizing symmetry and vector cancellation.

## How to Run

### Step 1: Install Python
- Ensure Python is installed on your system.

### Step 2: Install Dependencies
- Open a terminal or command prompt in this folder and run:
  ```bash
  pip install -r requirements.txt
  ```

### Step 3: Run the Simulation
- Copy the command below and run:
  python electric_field_viz.py

### Physics Orientation
- Charge Geometry: Vertical ring or disk in the Y–Z plane.
- Symmetry Axis: X-axis.
- Observer: Moves along the X-axis.

## Visual Features

### Vectors
- Yellow Arrows (dE): Individual electric field contributions.
- Scaling: VERY LARGE for visibility.
- Green Arrow (E_net): Net electric field at the observer.

### Toggles
- Located in the top-right corner.
- Moved to avoid overlap with sliders.

## Controls

### Sliders (Left Panel)
- X Position: Observer location.
- Radius (R): Ring or disk size.
- Charge (Q): Charge magnitude and sign.
- Count (N): Number of elements (vertical slider).

### Camera
- Orbit: Left mouse drag.
- Pan: Right mouse drag.
- Zoom: Scroll wheel.

### Camera Notes
- The camera orbits the target by default.
- The camera can be locked to the observer using "Follow Obs".

