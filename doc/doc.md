<style>
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
</style>

# MATLAB Code Documentation for Flow Field Post-processing

## Introduction
![alt text](<fig 1.png>)

 - A new step (Step 4) is introduced in this process to address the challenge of achieving convergence when only the maximum duration is considered. To overcome this limitation, we calculate the probability density function (PDF) of all durations during which the shear stress exceeds a predefined threshold. From this PDF, we derive a modified $t_m$ , which corresponds to a specific percentile (referred to as the "Given Ratio" in Figure 4 and as `data.ratio` in the code). Consequently, a revised threshold-duration relationship is established, as illustrated in Pic.5.
---

## Data Input Section
Below is a detailed explanation of the variables that need to be modified:

### Input Data File
- **`data_set`**: This variable specifies the name of the input data file which shall be in the same folder as the `post_processing.mat` file. The file should be in `.mat` format and contain the following variables:
  - `U`: Streamwise velocity component (3D array).
  - `V`: Spanwise velocity component (3D array).
  - `zpos_delta`: Wall-normal grid positions (1D array).
  - `xpos_delta`: Streamwise grid positions (1D array).

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
- Iterating over a set of ratios (`data.ratio`, which has been set ) to compute the critical value, physical duration, physical threshold and physical tau.

Plots of the critical droplet size for each ratio are generated and saved as PDF files in the `data_N` folder.

### Saving Results
The script saves the processed data (`data_N`) in a `.mat` file within the `data_N` folder. The file is named based on the Reynolds number.

## Output Files
After running the script, the following files will be generated in the `data_N` folder:
- **Mean Velocity Profile Plot**: A PDF file showing the mean velocity profile.
- **Critical Droplet Size Plots**: PDF files showing the critical droplet size for each ratio.
- **Processed Data File**: A `.mat` file containing the processed data.
