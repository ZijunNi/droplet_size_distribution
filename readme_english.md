<!-- <style>
  body {
    counter-reset: section; /* 初始化 section 计数器 */
  }
  h2 {
    counter-reset: subsection; /* 在每个 h2 中重置 subsection 计数器 */
  }
  h2::before {
    counter-increment: section; /* 递增 section 计数器 */
    content: counter(section) ". "; /* 显示 section 编号 */
  }
  h3::before {
    counter-increment: subsection; /* 递增 subsection 计数器 */
    content: counter(section) "." counter(subsection) ". "; /* 显示 subsection 编号 */
  }
</style> -->

# MATLAB Code Documentation for Flow Field Post-processing

若需要本文件的中文版本，请参见目录下的``readme.md``

## Introduction
![alt text](<./doc/fig 1.png>)

 - A new step (Step 4) is introduced in this process to address the challenge of achieving convergence when only the maximum duration is considered. To overcome this limitation, we calculate the probability density function (PDF) of all durations during which the shear stress exceeds a predefined threshold. From this PDF, we derive a modified $t_m$ , which corresponds to a specific percentile (referred to as the "Given Ratio" in Figure 4 and as `data.ratio` in the code). Consequently, a revised threshold-duration relationship is established, as illustrated in Pic.5.

## Data Input Section
Below is a detailed explanation of the variables that need to be modified:

### Input data File
- **`data_set`**: This variable specifies the name of the input data file which shall be in the same folder as the `post_processing.mat` file. The file should be in `.mat` format and contain the following variables:
  - `U`: Streamwise velocity component, normalized by the friction velocity $u_\tau$ (3D array).
  - `V`: Spanwise velocity component, normalized by the friction velocity $u_\tau$ (3D array).
  - `zpos_delta`: Wall-normal grid positions, normalized by the half width of the channel $\delta$ (1D array).
  - `xpos_delta`: Streamwise grid positions, normalized by the half width of the channel $\delta$ (1D array).

  Example:
  ```matlab
  data_set = "example_data.mat"; % Replace with the data file name
  ```

### Reynolds Number
- **`data.Reynolds_number`**: This variable specifies the Reynolds number of the flow. **This value must be consistent with the Reynolds number of the provided data**.

  Example:
  ```matlab
  data.Reynolds_number = 3200; % Replace with the Reynolds number
  ```
### The Upper and Lower Initial Values of the Bisection Method
- **`left_bound` and `right_bound`**: These two variables respectively represent the initial upper and lower bounds selected by the program when using the bisection method to approximate the critical value of particle diameter. **It is important to note that these two variables are normalized by the half-width of the channel.**


## Running the Script
Once the input data file and Reynolds number have been specified, you can run the script. The script will automatically perform the following steps:

### Loading Data
The script loads the input data file and extracts the necessary variables (`U`, `V`, `zpos_delta`, `xpos_delta`).

### Creating Data Folder
The script creates a folder named `data_N` in the current working directory to store the results and plots, where `N` stands for the certain Reynolds number `data.Reynolds_number`. If the folder already exists, it will not be recreated.

### Calculating Mean Velocity Profile
The script calculates the mean velocity profile by averaging the streamwise velocity component (`U`) across the spanwise and streamwise directions. The result is saved in `data.mean_U`.

A plot of the mean velocity profile is generated and saved as a PDF file in the `data_N` folder.

### Calculating Critical Droplet Size
The script calculates the critical droplet size based on the input data. The calculation involves the following steps:
- Defining the bounds for the initial droplet size (`left_bound` and `right_bound`).
- Iterating over a set of ratios (`data.ratio`, which has been set ) to compute the critical value, physical duration, physical threshold and physical tau of droplet.

Plots of the critical droplet size for each ratio are generated and saved as PDF files in the `data_N` folder.

### Saving Results
The script saves the processed data (`data_N`) in a `.mat` file within the `data_N` folder. The file is named based on the Reynolds number.

## Output Files
After running the script, the following files will be generated in the `data_N` folder:
- **Mean Velocity Profile Plot**: A PDF file showing the mean velocity profile.
- **Critical Droplet Size Plots**: PDF files showing the critical droplet size for each ratio.
- **Processed Data File**: A `.mat` file containing the processed data.

## Error and Warning Notifications
- **<font color=red>Initial Bound Error, [left_bound, right_bound] need to be adjusted.</font>**
This notification suggests that the current upper and lower bounds are not appropriate for the dataset being used. The issue can be resolved by adjusting the two variables described in **Section 2.3**. As a practical reference, the critical value of particle diameter is generally slightly smaller than $y^+\approx 1$, which can guide the adjustment process.
If the aforementioned method does not yield the desired results, consider modifying the parameter within the parentheses following the `if` statement on line 8 of the file `extract_2d_slice_x_interp.m`. By default, this parameter is set to 1, which enables logarithmic interpolation to obtain the results rather than simple linear interpolation. However, in cases where both low friction Reynolds numbers and wall-resolved grids are present, the near-wall grid points may already lie within the linear region. In such scenarios, linear interpolation should be used instead, which can be achieved by setting the aforementioned parameter to 0. Please note that if this value is modified, it is essential to revert it back to the default value of 1 before running the program again.

- **<font color=red>The mean velocity profile deviates from the logarithmic law. Please verify that the input velocity values are properly non-dimensionalized.</font>**
The aforementioned notification indicates that the mean velocity profile of the input data deviates significantly from the standard logarithmic law. **It is important to note that the input velocity field should be non-dimensionalized using the friction velocity $u_\tau$, while the input coordinates should be normalized by the half-width of the channel $\delta$.** The logarithmic law referenced here is given by $u^+ = \ln(y^+)/\kappa + B$, where $\kappa = 0.41$ and $B = 5$.

## 6. Appendix: Detailed Function Specifications

### 6.1 `critical_value` Function
**Input Parameters Table**  
| Parameter Name | Parameter Meaning | Non-dimensional Units |
| :--: | :--: | :--: |
| `U`, `V` | Streamwise and spanwise components of velocity field | Friction velocity $u_\tau$ |
| `dx` | Streamwise grid spacing | Half-channel width $\delta$ |
| `zpos_delta` | Vertical grid position | Half-channel width $\delta$ |
| `Reynolds_number` | Friction velocity-based Reynolds number | \ |
| `ratio` | Quantile threshold ratio | \ |
| `left_bound`, `right_bound` | Initial bounds for bisection method | \ |

The function uses a constant continuous phase density of 1000 and employs `condition_function` to determine whether droplets of specified size will break.

### 6.2 `condition_function` Function

**Input Parameters Table**  
| Parameter Name | Parameter Meaning | Non-dimensional Units |
| :--: | :--: | :--: |
| `u`, `v` | Streamwise and spanwise components of velocity field | Friction velocity $u_\tau$ |
| `dx` | Streamwise grid spacing | Half-channel width $\delta$ |
| `zpos_delta` | Vertical grid position | Half-channel width $\delta$ |
| `diameter` | Particle diameter | Half-channel width $\delta$ |
| `Reynolds_number` | Friction velocity-based Reynolds number | \ |
| `ratio` | Quantile threshold ratio | \ |
| `rho_c` | Continuous phase density | With units |

**Output Parameters Table**  
| Parameter Name | Type | Description |
|----------------|------|-------------|
| `result` | Boolean | Whether erosion condition is met (True/False) |
| `physical_duration` | Array | Effective stress duration (unit: seconds) |
| `physical_threshold` | Array | Physical critical shear stress threshold (unit: Pa) |
| `physical_tau` | Array | Particle breakup characteristic time (unit: seconds) |

#### Algorithm Workflow
1. **Physical Quantity Calculation**
   - Compute friction velocity using experimental parameters (half-channel width $\delta = 5\ \rm{mm}$ and continuous phase viscosity $\nu = 1.8\times 10^{-6}\ \rm{m^2/s}$) and input Reynolds number:  
     $u_\tau = Re \cdot \nu_c / \delta$
   - Calculate non-dimensional stress threshold `threshold_series` and corresponding duration `duration` using `threshold_line` function.

2. **Unit Conversion**  
   - Convert non-dimensional stress to physical units using $\tau = \rho_c u_d^2 / 2$:  
     `physical_threshold = threshold_series * 0.5 * rho_c * u_tau^2`
   - Convert spatial scales to temporal scales via Taylor's frozen hypothesis. Three velocity estimation methods are implemented:  
     - Local mean velocity  
     - Fixed multiple (0.8×) of centerline velocity  
     - Fixed multiple (9.5×) of friction velocity  
     Compute $\bar{u}$ and convert duration:  
     `physical_duration = duration * delta / u_bar`

3. **Breakup Condition Judgment**  
   - Compute particle breakup characteristic time `physical_tau` using `calculate_tau` function.
   - Determine breakup condition via `physical_tau < physical_threshold` and output result in `result`.

### 6.3 `calculate_tau` Function

**Input Parameters Table**  
| Parameter Name | Parameter Meaning | Non-dimensionalized Unit |  
| :--: | :--: | :--: |  
| `t_c` | Duration of shear stress acting on the droplet | Dimensional (with units) |  
| `d` | Droplet diameter | Half-channel width δ |  

**Predefined Physical Parameters in the Function**  
| Parameter Name | Parameter Meaning | Value |  
| :--: | :--: | :--: |  
| `mu_d` | Droplet fluid viscosity | 0.0021 kg/(m·s) |  
| `sigma` | Surface tension coefficient | 0.004 N/m |  
| `theta_c` | Critical deformation threshold | 1 |  

**Implementation Notes**  
The function provides two methods to calculate the critical breakup criteria for droplets. Generally, the nonlinear term can be neglected, and the first method suffices. The second model directly derives from the original Voigt model, requiring symbolic solutions to differential equations and is more complex.  

**Function Description**  
For each input threshold duration, the function calculates and outputs the maximum shear stress `tau` that the droplet can withstand without breaking under that duration.