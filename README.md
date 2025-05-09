# Thermal Conductivity Calculation

A numerical method for predicting the effective thermal conductivity of filler-filled polymer composites. This approach considers:
- Filler properties (thermal conductivity, particle radius)
- Polymer matrix properties (thermal conductivity)
- Composite characteristics (filler volume fraction, interfacial thermal resistance)

## System Requirements

### Python Environment
- **Minimum**: Python 3.9+

### Required Packages
```bash
pip install flask sympy joblib scipy numpy matplotlib pandas
```

## Installation
1. Download the repository:
   - Click the green "Code" button
   - Select "Download ZIP"
   - Save as `Thermal-conductivity-calculation-main.zip`

Alternative direct download:
```bash
wget https://github.com/ThermalCalc/Thermal-conductivity-calculation/archive/refs/heads/main.zip
```

## Usage

### Web Interface
1. Run the server:
   ```bash
   python server.py
   ```
2. Access the interface at: `http://127.0.0.1:5000`

### Calculation Modules

#### 1. Single Dataset Calculation
Enter parameters directly in the web form:
- Filler properties
- Polymer properties
- Composite parameters

#### 2. Batch Processing
1. Prepare input file (`data.txt`):
   ```
   volume_fraction, particle_size, polymer_conductivity, filler_conductivity, interface_resistance
   0.1, 1e-6, 0.2, 100, 1e-8
   0.2, 1e-6, 0.2, 100, 1e-8
   ```
2. Upload via "Batch Data Set Calculation" module
3. Results saved in `/results` directory

#### 3. Interface Resistance Calculation
Select model (EMT/BCC) and input:
- Composite thermal conductivity
- Constituent material properties
- Filler characteristics

## MATLAB Implementation
The repository includes:
- `Effective_Thermal_Conductivity_Calculation.m` script
- Separate implementations for:
  - Low-fill systems (Vf ≤ 13%)
  - High-fill systems (13% < Vf ≤ 68%)

### MATLAB Usage
1. Prepare input file (`xx.data`) formatted as specified in lines 5-13
2. Alternatively, modify parameters directly in the script
3. Run to obtain:
   - Predicted effective thermal conductivity
   - Model validation against experimental data (if available)


For detailed theoretical background, please refer to our [documentation](docs/theory.pdf).
```
