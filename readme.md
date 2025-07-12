# 2D Ray-Tracing Simulator for V2V Channel Modeling

## 1. Project Overview

This MATLAB project provides a 2D ray-tracing simulator designed to model wireless signal propagation in environments with obstacles, such as Vehicle-to-Vehicle (V2V) communication scenarios. The core of the simulator is the **Method of Images**, a powerful algorithm used to find all possible propagation paths—including the direct Line-of-Sight (LOS) path and multiple reflected paths—between a transmitter (TX) and a receiver (RX).

For each valid path found, the simulator calculates the **complex channel gain (α)**. This value is critical in wireless communications as it encapsulates both the amplitude attenuation (path loss and reflection losses) and the phase shift the signal undergoes. The collection of all path gains forms the channel impulse response, which is fundamental for analyzing wireless channel performance.

The project is structured in a modular way, with each function having a single, well-defined responsibility. This makes the code clean, easy to understand, and scalable.

## 2. How to Run the Simulation

1.  **File Placement**: Ensure all the `.m` function files are located in the same directory as the main script, or that the folder containing them is added to the MATLAB path.
2.  **Configure the Scenario**: Open the main script `main.m`. Inside this file, you can easily modify the simulation scenario:
    * Define the environment by changing the coordinates of the `walls`.
    * Set the `transmitter_position` and `receiver_position`.
    * Adjust `k_max` to control the maximum number of reflections to compute.
    * Modify physical constants in the `sim_params` struct.
3.  **Execute**: Run the `main.m` script from the MATLAB command window or editor.
4.  **View Results**: A plot will be generated showing the environment, the TX and RX, and all the valid propagation paths found by the simulator. The MATLAB console will print the total number of paths found.

## 3. Core Methodology: The Method of Images

The algorithm is built on the **Method of Images**, a geometric technique to find reflected paths. The core idea is that a reflected ray can be seen as a straight-line ray originating from a "virtual" source, which is the mirror image of the original source.

* **1st-Order Reflection**: To find a path that reflects off one wall, we first find the image of the transmitter, `TX'`, by reflecting its position across the line defined by the wall. A straight line is then drawn from this image `TX'` to the receiver `RX`. The point where this line intersects the wall is the reflection point `R1`. The valid path is `TX -> R1 -> RX`.
* **Higher-Order Reflections**: The process is extended recursively. To find a 2nd-order path (reflecting off Wall 1, then Wall 2), we first find the image `TX'` across Wall 1. Then, we reflect `TX'` across Wall 2 to get a new image, `TX''`. The path is found by tracing back from `RX` to `TX''`.

This recursive creation of images allows the algorithm to systematically find all possible reflection sequences up to a given order `k_max`.

## 4. Function by Function Breakdown

The code is broken down into several specialized functions. Here is a detailed look at each one.

### `main.m`

* **Purpose**: The main entry point of the simulation. It is a script, not a function.
* **Logic**:
    1.  Clears the workspace (`clear`), command window (`clc`), and closes all figures (`close all`).
    2.  Sets up the simulation `sim_params` struct, which holds all physical constants.
    3.  Defines the `walls` of the environment and the `tx_pos` and `rx_pos`.
    4.  Calls the main computational engine, `runRayTracing`, to get all path data.
    5.  Passes the results to `plotRays` for visualization.

### `runRayTracing.m`

* **Purpose**: The high-level manager of the calculation process.
* **Logic**:
    1.  **LOS Path**: It first checks for the existence of a direct Line-of-Sight path. It does this by calling `findSegmentIntersection` between the TX-RX line and every wall. If no intersections are found, a valid LOS path exists, and its properties are calculated and stored.
    2.  **Reflected Paths**: It then loops from `order = 1` to `k_max`. In each iteration, it calls the core recursive function `findReflectedraysRecursive` to find all paths with that specific number of reflections.
    3.  **Data Aggregation**: It collects the results from the LOS check and all reflection orders into a single cell array, `all_rays_data`.
    4.  **Output Formatting**: It extracts just the complex gains (`alphas`) into a simple vector for convenience.

### `findReflectedraysRecursive.m`

* **Purpose**: The heart of the image method algorithm. It recursively explores the "tree" of all possible reflections.
* **Logic**:
    * **Base Case (`reflections_remaining == 0`)**: This is the termination condition for the recursion. When a full sequence of images has been generated, this function calls `validaterayGeometry`. If the geometry is valid (i.e., a real physical path exists), it then calls `calculatePhysicalProperties` to compute the path's details and returns the resulting data structure.
    * **Recursive Step**: If more reflections are needed, the function iterates through all `walls`. For each wall, it creates a new image source by calling `reflectPointAcrossWall`. It then calls **itself** with this new image source, decrementing `reflections_remaining` and updating the history of walls and images used. This branching process builds the entire tree of possibilities.

### `validaterayGeometry.m`

* **Purpose**: To verify if a theoretical path from a sequence of images corresponds to a real, unobstructed physical path.
* **Logic**:
    1.  **Trace Backwards**: The function starts from the `receiver_position` and traces back towards the final image in the sequence.
    2.  **Find Reflection Point**: The intersection of the line `(Image -> RX)` and the corresponding wall gives the first reflection point. It calls `findSegmentIntersection` for this.
    3.  **Validation 1 (On Segment)**: It checks if this intersection point lies on the finite physical wall segment. If not, the path is invalid.
    4.  **Validation 2 (Obstruction)**: It then checks if the newly traced segment `(Reflection Point -> RX)` is blocked by any *other* walls. If it is, the path is invalid.
    5.  **Iterate**: The process repeats, tracing back from the newly found reflection point to the next image in the sequence, until the full path to the original `transmitter_position` is traced and validated. If all checks pass, it returns the coordinates of the path vertices.

### `calculatePhysicalProperties.m`

* **Purpose**: To compute the physical characteristics of a geometrically valid path.
* **Logic**:
    1.  **Distance**: It calculates the total path length by summing the Euclidean distance of each segment.
    2.  **Reflection Loss**: It iterates through each reflection point. For each reflection, it calculates the **Fresnel Reflection Coefficient (Γ)**. This complex number represents the ratio of reflected to incident electric field amplitude. The formula used is for perpendicular polarization:
        $$ 
        \Gamma_{\perp} = \frac{\cos(\theta_i) - \sqrt{\varepsilon_r - \sin^2(\theta_i)}}{\cos(\theta_i) + \sqrt{\varepsilon_r - \sin^2(\theta_i)}} 
        $$
        where $ \theta_i $ is the angle of incidence and $ \varepsilon_r $ is the relative permittivity of the wall material.
    3.  **Cumulative Product**: It multiplies all the individual reflection coefficients together to get `cumulative_gamma` ($ \prod \Gamma_i $).
    4.  **Final Gain**: It calls `calculateAlpha_n` to compute the final complex gain using the total distance and cumulative reflection loss.

### `calculateAlpha_n.m`

* **Purpose**: To calculate the final complex channel gain `alpha_n` for a single path.
* **Logic**: This function implements the channel gain formula, derived from the Friis transmission equation and wave propagation principles.
    $$ \alpha_n = j \frac{\lambda Z_0}{4\pi^2 R_a d_n} \left( \prod_{i=1}^{N} \Gamma_i \right) e^{-j 2\pi f_c \tau_n} $$
    * **Amplitude Term**: $ \frac{\lambda Z_0}{4\pi^2 R_a d_n} $ accounts for free-space path loss (power spreading over distance $d_n$) and antenna/impedance properties.
    * **Reflection Term**: $ \prod_{i=1}^{N} \Gamma_i $ is the `cumulative_gamma` calculated previously, representing the total loss from all reflections.
    * **Phase Term**: $ e^{-j 2\pi f_c \tau_n} $ represents the phase shift of the signal due to the travel time (delay $ \tau_n = d_n / c $) along the path.

### `reflectPointAcrossWall.m`

* **Purpose**: A geometric utility to reflect a point across a line.
* **Logic**: It uses the standard vector formula for reflection. For a point `P` to be reflected across a line that passes through point `L` and has a normal vector `n`, the reflected point `P'` is:
    $$ P' = P - 2 \cdot \frac{(P - L) \cdot \mathbf{n}}{||\mathbf{n}||^2} \cdot \mathbf{n} $$

### `findSegmentIntersection.m`

* **Purpose**: A geometric utility to find the intersection point of two line segments.
* **Logic**: It uses a vector cross-product method. Given two segments `P1 = P1_a + t(P1_b - P1_a)` and `P2 = P2_a + u(P2_b - P2_a)`, it solves for the parameters `t` and `u`. An intersection exists only if both `t` and `u` are between 0 and 1.

### `plotRays.m`

* **Purpose**: To create a visual representation of the simulation results.
* **Logic**:
    1.  Creates a new figure and sets up the axes for an equal-aspect-ratio plot.
    2.  Plots the `walls` as black lines.
    3.  Plots the `transmitter_position` and `receiver_position` as colored circles.
    4.  Loops through the `all_rays_data`. For each path, it plots the series of vertices.
    5.  It uses a different color for each reflection order (e.g., green for LOS, blue for 1-reflection, red for 2-reflections, etc.) to make the plot easy to interpret.
    6.  Finally, it adds a title, axis labels, and a clean legend.