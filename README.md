# Thermal-conductivity-calculation
Calculation of effective thermal conductivity of complex compositions.

The **Thermal Conductivity Calculation** method provides a numerical approach to accurately predict the effective thermal conductivity of polymer composites. This method considers various parameters, including the properties of the filler (thermal conductivity and particle size), the properties of the polymer (thermal conductivity), and the properties of the composite material (volume fraction of the filler and interfacial thermal resistance between the filler and the polymer).

## Steps to Operate:

1. **Run the PYTHON Software**:
   - Open the PYTHON software and follow the instructions to open the web interface.

2. **Input Parameters**:
   - Enter the required parameters directly into the provided fields and click the button to calculate the effective thermal conductivity.

3. **Batch Calculations for Multiple Data Sets**:
   - For calculating multiple sets of data, create a `.txt` file with the following format:
     - **First Column**: Volume fraction of the filler
     - **Second Column**: Particle size of the filler
     - **Third Column**: Thermal conductivity of the polymer
     - **Fourth Column**: Thermal conductivity of the filler
     - **Fifth Column**: Interfacial thermal resistance between the filler and the polymer

4. **Interfacial Thermal Resistance Calculation**:
   - The interfacial thermal resistance between the filler and the polymer can be calculated using either the EMT model or the BC model, depending on the requirements.
