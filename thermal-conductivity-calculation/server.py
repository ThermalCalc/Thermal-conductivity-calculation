from flask import Flask, render_template, request, jsonify, send_from_directory
from sympy import symbols, sqrt, pi, integrate, lambdify,simplify,hyper,Abs,re,im
from joblib import Parallel, delayed
from scipy.integrate import quad
from scipy.special import hyp2f1
from math import isfinite

import os
import json
import numpy as np
import pandas as pd


# Initialize Flask Application
app = Flask(__name__, static_folder='static')  
UPLOAD_FOLDER = './uploads'
RESULT_FOLDER = './results'

# Create Necessary Folders
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(RESULT_FOLDER, exist_ok=True)


# Define Symbolic Variables**
x, m = symbols('x m', real=True)

# Function to Calculate Keff for a Single Data Set**
def calculate_keff_single(r, Vf, Kf, Km, Ri):
    """
    Calculate Keff for a single data set based on the formula.
    """
    try:
        if Vf<0.131:
            Keff=calculate_keff_low_Vf(r, Vf, Kf, Km, Ri)
        else:
            Keff=calculate_keff_high_Vf(r, Vf, Kf, Km, Ri)
        return {"Keff": Keff}
    except Exception as e:
        print(f"Error during calculation: {e}")
        return {"error": str(e)}

def calculate_keff_low_Vf(r, Vf, Kf, Km, Ri):  
    try:       
        l = (r * (pi / 3 / Vf) ** (1 / 3)).evalf()
        a = (sqrt(2) * r).evalf()
        b = (sqrt(2 / 3) * r).evalf()
        w = 2*r/sqrt(3).evalf()

        # Calculate thermal resistance R1.
        ex_expr=(-3*m+sqrt(-3*m**2+4*r**2))/(2*sqrt(2))
        ey_expr=(m+sqrt(-3*m**2+4*r**2))/(2*sqrt(2))
        #0<m<r
        def A1_R1(m_val):
            m_val = float(m_val)
            ex_value = float(ex_expr.subs(m, m_val).evalf())
            ey_value = float(ey_expr.subs(m, m_val).evalf())
            def integrand(x):
                return b / a * (a ** 2 - x ** 2) ** 0.5
            integral_value1, _ = quad(integrand, ex_value, a)
            return ey_value**2+2*integral_value1    
        def integrand_R1(m_val):
            try:
                A1_val = float(A1_R1(m_val)) 
                Vf1 = float(pi * (r**2 - m_val**2) / (4 * A1_val))
                Vm1 = 1 - Vf1
                K1 = float(Vf1 * Kf + Vm1 * Km)
                result_1 = float(1 / (A1_val * K1))

                return result_1

            except Exception as e:
                print(f"Error in integrand_R1 for m_val={m_val}: {e}")
                raise  # Re-raise the exception for debugging purposes
        
        # Calculate R1 using scipy.integrate.quad.
        R1, _ = quad(integrand_R1, 0, r)

        # Calculate thermal resistance R2.
        kx_expr=(-3*m-sqrt(-3*m**2+4*r**2))/(2*sqrt(2))
        ky_expr=(m-sqrt(-3*m**2+4*r**2))/(2*sqrt(2))
        #r<m<w
        def A2_R2(m_val):
            ex_value = float(ex_expr.subs(m, m_val).evalf())
            kx_value = float(kx_expr.subs(m, m_val).evalf())
            ey_value = float(ey_expr.subs(m, m_val).evalf())
            ky_value = float(ky_expr.subs(m, m_val).evalf())
            
            def integrand(x):
                return b / a * (a ** 2 - x ** 2) ** 0.5
            integral_value2, _ = quad(integrand, -a, kx_value)
            integral_value3, _ = quad(integrand, ex_value, a)
            return 2 * integral_value2+(ey_value+ky_value)*(ex_value-kx_value)+2*integral_value3

        # Define the numerical integration of R2
        def integrand_R2_1(m_val):
            try:
                A2_val = float(A2_R2(m_val))
                K2 = float(Km)
                result_2 = float(1 / (A2_val * K2))
                return result_2

            except Exception as e:
                raise
            
        R2_1, _ = quad(integrand_R2_1, r, w)
        #w<m<l-w
        R2_2=(l-2*w)/(pi*a*b*Km)
        
        R2=R2_1+R2_2
        
        # Calculate thermal resistance R3.
        R3=R1
        
        # Calculate thermal resistance R4.
        #0<m<r
        def integrand_R4_1(m_val):
            try:
                A4_val1 = float(l**2-A1_R1(m_val))
                result_R4_1 = float(1 / (A4_val1 * Km))
                return result_R4_1
            except Exception as e:
                print(f"Error in integrand_R4_1 for m_val={m_val}: {e}")
                raise
        R4_1, _ = quad(integrand_R4_1, 0, r)
        #r<m<w
        def integrand_R4_2(m_val):
            try:
                A4_val2 = float(l**2-A2_R2(m_val))
                result_R4_2 = float(1 / (A4_val2 * Km))
                return result_R4_2

            except Exception as e:
                print(f"Error in integrand_R4_1 for m_val={m_val}: {e}")
                raise
        R4_2, _ = quad(integrand_R4_2, r, w)
        #w<m<l-w
        R4_3 = (l-2*w)/((l**2-pi*a*b)*Km)
        R4 = 2*R4_1+2*R4_2+R4_3

        # Calculate the interfacial thermal resistance, RI.
        RI = float(Ri / (pi * r ** 2 / 2))

        # Calculate the thermal resistance of the minimum structural unit, Re.
        Re = float(1 / (1 / (2 * R1 + R2 + 2 * RI) + 1 / R4))

        # Calculate the effective thermal conductivity, Keff.
        Keff = float( l / l ** 2 /2/ Re)
        return {"Keff": Keff}
    except Exception as e:
        print(f"Error during calculation: {e}")
        return {"error": str(e)}

def calculate_keff_high_Vf(r, Vf, Kf, Km, Ri):  
    try:       
        l = (r * (pi / 3 / Vf) ** (1 / 3)).evalf()
        a = (sqrt(2) * r).evalf()
        b = (sqrt(2 / 3) * r).evalf()
        w = 2*r/sqrt(3).evalf()

        # Calculate thermal resistance R1.
        ex_expr=(-3*m+sqrt(-3*m**2+4*r**2))/(2*sqrt(2))
        ey_expr=(m+sqrt(-3*m**2+4*r**2))/(2*sqrt(2))
        mx_expr=(3*l-3*m-sqrt(-3*(l-m)**2+4*r**2))/(2*sqrt(2))
        my_expr=(l-m+sqrt(-3*(l-m)**2+4*r**2))/(2*sqrt(2))
        ox_expr=(3*l-3*m+sqrt(-3*(l-m)**2+4*r**2))/(2*sqrt(2))
        oy_expr=(l-m-sqrt(-3*(l-m)**2+4*r**2))/(2*sqrt(2))
        #0<m<l-w
        def A1_R1_1(m_val):
            m_val = float(m_val)
            ex_value = float(ex_expr.subs(m, m_val).evalf())
            ey_value = float(ey_expr.subs(m, m_val).evalf())
            def integrand(x):
                return b / a * (a ** 2 - x ** 2) ** 0.5
            integral_value1, _ = quad(integrand, ex_value, a)
            return ey_value**2+2*integral_value1
        def integrand_R1_1(m_val):
            try:
                A1_val1 = float(A1_R1_1(m_val)) 
                Vf1 = float(pi * (r**2 - m_val**2) / (4 * A1_val1))
                Vm1 = 1 - Vf1
                K1 = Vf1 * Kf + Vm1 * Km
                K1 = float(K1)
                result_1 = 1 / (A1_val1 * K1)
                result_1 = float(result_1)
                return result_1

            except Exception as e:
                print(f"Error in integrand_R1 for m_val={m_val}: {e}")
                raise  # Re-raise the exception for debugging purposes
        #l-w<m<l-r    
        def A1_R1_2(m_val):
            m_val = float(m_val)
            ex_value = float(ex_expr.subs(m, m_val).evalf())
            ey_value = float(ey_expr.subs(m, m_val).evalf())
            ox_value = float(ox_expr.subs(m, m_val).evalf())
            oy_value = float(oy_expr.subs(m, m_val).evalf())
            mx_value = float(mx_expr.subs(m, m_val).evalf())
            my_value = float(my_expr.subs(m, m_val).evalf())
            def integrand(x):
                return b / a * (a ** 2 - x ** 2) ** 0.5
            integral_value2, _ = quad(integrand, ex_value, mx_value)
            integral_value3, _ = quad(integrand, ox_value, a)
            return ey_value**2+2*integral_value2+(oy_value+my_value)*(ox_value-mx_value)+2*integral_value3
        def integrand_R1_2(m_val):
            try:
                A1_val2 = float(A1_R1_2(m_val)) 
                Vf1 = float(pi * (r**2 - m_val**2) / (4 * A1_val2))
                Vm1 = 1 - Vf1
                K1 = float(Vf1 * Kf + Vm1 * Km)
                result_2 = float(1 / (A1_val2 * K1))
                return result_2
            except Exception as e:
                print(f"Error in integrand_R1 for m_val={m_val}: {e}")
                raise  # Re-raise the exception for debugging purposes
        
        # Calculate R1 using scipy.integrate.quad.
        R1_1, _ = quad(integrand_R1_1, 0, l-w)
        R1_2, _ = quad(integrand_R1_2, l-w, l-r)
        R1 = R1_1+R1_2

        # Calculate thermal resistance R2.
        #l-r<m<r
        def A2_R2(m_val):
            ex_value = float(ex_expr.subs(m, m_val).evalf())
            mx_value = float(mx_expr.subs(m, m_val).evalf())
            ey_value = float(ey_expr.subs(m, m_val).evalf())
            my_value = float(my_expr.subs(m, m_val).evalf())
            def integrand(x):
                return b / a * (a ** 2 - x ** 2) ** 0.5
            integral_value, _ = quad(integrand, ex_value, mx_value)
            return ey_value**2+2 * integral_value+my_value**2
        def integrand_R2(m_val):
            try:
                A2_val = float(A2_R2(m_val))
                Vf2 = float(pi * (r**2 - m_val**2+r**2-(l-m_val)**2) / (4 * A2_val))
                Vm2 = 1 - Vf2
                K2 = float(Vf2 * Kf + Vm2 * Km)
                result_2 = float(1 / (A2_val * K2))
                return result_2
            except Exception as e:
                raise
        R2, _ = quad(integrand_R2, l-r, r)
        
        # Calculate thermal resistance R4.
        #0<m<l-w
        def integrand_R4_1(m_val):
            try:
                A4_val1 = float(l**2-A1_R1_1(m_val))
                result_R4_1 = float(1 / (A4_val1 * Km))
                return result_R4_1
            except Exception as e:
                print(f"Error in integrand_R4_1 for m_val={m_val}: {e}")
                raise
        R4_1, _ = quad(integrand_R4_1, 0, l-w)
        #l-w<m<l-r
        def integrand_R4_2(m_val):
            try:
                A4_val2 = float(l**2-A1_R1_2(m_val))
                result_R4_2 = float(1 / (A4_val2 * Km))
                return result_R4_2
            except Exception as e:
                print(f"Error in integrand_R4_1 for m_val={m_val}: {e}")
                raise
        R4_2, _ = quad(integrand_R4_2, l-w, l-r)
        #l-r<m<r
        def integrand_R4_3(m_val):
            try:
                A4_val3 = float(l**2-A2_R2(m_val))
                result_R4_3 = float(1 / (A4_val3 * Km))
                return result_R4_3
            except Exception as e:
                print(f"Error in integrand_R4_1 for m_val={m_val}: {e}")
                raise
        R4_3,_ = quad(integrand_R4_2, l-r, r)
        R4 = 2*R4_1+2*R4_2+R4_3

        # Calculate the interfacial thermal resistance, RI.
        RI = float(Ri / (pi * r ** 2 / 2))

        # Calculate the thermal resistance of the minimum structural unit, Re.
        Re = float(1 / (1 / (2 * R1 + R2 + 2 * RI) + 1 / R4))

        # Calculate the effective thermal conductivity, Keff.
        Keff = float( l / l ** 2 /2/ Re)
        
        return {"Keff": Keff}

    except Exception as e:
        print(f"Error during calculation: {e}")
        return {"error": str(e)}

# Route: Home page
@app.route('/')
def index():
    return send_from_directory(app.static_folder, 'Web-based calculation design.html')

# Route: Single Data Set Calculation
@app.route('/calculate-single', methods=['POST'])
def calculate_single():
    try:
        data = request.json

        r = float(data['r'])
        Vf = float(data['Vf'])
        Kf = float(data['Kf'])
        Km = float(data['Km'])
        Ri = float(data['Ri'])

        result = calculate_keff_single(r, Vf, Kf, Km, Ri)
        return jsonify(result)
    except Exception as e:
        return jsonify({"error": str(e)})

@app.route('/calculate-multiple', methods=['POST'])
def calculate_multiple():
    try:
        # Retrieve the uploaded file
        file = request.files['file']
        file_path = os.path.join(UPLOAD_FOLDER, file.filename)
        file.save(file_path)

        # Load the file data (assuming the file is .txt or .csv, with data separated by commas)
        data = np.loadtxt(file_path, delimiter=',')  # Data format check

        # Data columns: r, Vf, Kf, Km, Ri
        r = data[:, 0]
        Vf = data[:, 1]
        Kf = data[:, 2]
        Km = data[:, 3]
        Ri = data[:, 4]

        # Parallel computation
        results = Parallel(n_jobs=-1)(
            delayed(calculate_keff_single)(r[i], Vf[i], Kf[i], Km[i], Ri[i]) for i in range(len(r))
        )

        # Save the calculation results
        result_path = os.path.join(RESULT_FOLDER, 'results.txt')
        with open(result_path, 'w') as f:
            for result in results:
                if "error" in result:
                    f.write(f"Error: {result['error']}\n")
                else:
                    f.write(f"{result['Keff']}\n")

        # Return a JSON response with the download path.
        return jsonify({
            "success": True,
            "result": {
                "filePath": result_path,
                "downloadUrl": f"/download-results/{os.path.basename(result_path)}"
            }
        })
    except Exception as e:
        # Catch exceptions and return error messages.
        return jsonify({"success": False, "error": str(e)}), 500


@app.route('/download-results/<filename>', methods=['GET'])
def download_results(filename):
    try:
        return send_from_directory(RESULT_FOLDER, filename, as_attachment=True)
    except FileNotFoundError:
        return jsonify({"success": False, "error": "File not found"}), 404


# Calculate interface resistance
@app.route('/interface-resistance')
def interface_page():
    return render_template('interface.html')

# Route for interface resistance calculation
@app.route('/calculate-interface', methods=['POST'])
def calculate_interface():
    try:
        data = request.json
        Ka = float(data['Ka'])
        Vf = float(data['Vf'])
        Kf = float(data['Kf'])
        Km = float(data['Km'])
        r = float(data['r'])
        model = data['model']

        if model == 'emt':
            
            # Example EMT calculation
            Ri = r / Km *(( Kf - Km)/ Kf+3*(Km - Ka)/(2 * Km * (Vf - 1)+ Ka * (2 + Vf)))
            
        elif model == 'bc':

            # Example BC calculation
            m = symbols('m', real=True)
            l = (r * (pi / 3 / Vf) ** (1 / 3)).evalf()
            a = (sqrt(2) * r).evalf()
            b = (sqrt(2 / 3) * r).evalf()
            w = 2*r/sqrt(3).evalf()
            if Vf<0.131:
                # Calculate thermal resistance R1.
                ex_expr=(-3*m+sqrt(-3*m**2+4*r**2))/(2*sqrt(2))
                ey_expr=(m+sqrt(-3*m**2+4*r**2))/(2*sqrt(2))
    
                def A1_R1(m_val):
                    m_val = float(m_val)
                    ex_value = float(ex_expr.subs(m, m_val).evalf())
                    ey_value = float(ey_expr.subs(m, m_val).evalf())
                    def integrand(x):
                        return b / a * (a ** 2 - x ** 2) ** 0.5
                    integral_value1, _ = quad(integrand, ex_value, a)
                    return ey_value**2+2*integral_value1
                def integrand_R1(m_val):
                    try:
                        A1_val = float(A1_R1(m_val))
                        Vf1 = float(pi * (r**2 - m_val**2) / (4 * A1_val))
                        Vm1 = 1 - Vf1
                        K1 = float(Vf1 * Kf + Vm1 * Km)
                        result_1 = float(1 / (A1_val * K1))
                        return result_1
                    except Exception as e:
                        print(f"Error in integrand_R1 for m_val={m_val}: {e}")
                        raise  # Re-raise the exception for debugging purposes
                R1, _ = quad(integrand_R1, 0, r)

                # Calculate thermal resistance R2.
                kx_expr=(-3*m-sqrt(-3*m**2+4*r**2))/(2*sqrt(2))
                ky_expr=(m-sqrt(-3*m**2+4*r**2))/(2*sqrt(2))
                def A2_R2(m_val):
                    ex_value = float(ex_expr.subs(m, m_val).evalf())
                    kx_value = float(kx_expr.subs(m, m_val).evalf())
                    ey_value = float(ey_expr.subs(m, m_val).evalf())
                    ky_value = float(ky_expr.subs(m, m_val).evalf())
                    def integrand(x):
                        return b / a * (a ** 2 - x ** 2) ** 0.5
                    integral_value2, _ = quad(integrand, -a, kx_value)
                    integral_value3, _ = quad(integrand, ex_value, a)
                    return 2 * integral_value2+(ey_value+ky_value)*(ex_value-kx_value)+2*integral_value3
                def integrand_R2_1(m_val):
                    try:
                        A2_val = float(A2_R2(m_val))
                        K2 = float(Km)
                        result_2 = float(1 / (A2_val * K2))
                        return result_2
                    except Exception as e:
                        raise
                R2_1, _ = quad(integrand_R2_1, r, w)
                R2=R2_1+(l-2*w)/(pi*a*b*Km)
                # Calculate thermal resistance R4.
                def integrand_R4_1(m_val):
                    try:
                        A4_val1 = float(l**2-A1_R1(m_val))
                        result_R4_1 = float(1 / (A4_val1 * Km))
                        return result_R4_1
                    except Exception as e:
                        print(f"Error in integrand_R4_1 for m_val={m_val}: {e}")
                        raise
                R4_1, _ = quad(integrand_R4_1, 0, r)
                def integrand_R4_2(m_val):
                    try:
                        A4_val2 = float(l**2-A2_R2(m_val))
                        result_R4_2 = float(1 / (A4_val2 * Km))
                        return result_R4_2
                    except Exception as e:
                        print(f"Error in integrand_R4_1 for m_val={m_val}: {e}")
                        raise
                R4_2, _ = quad(integrand_R4_2, r, w)
                R4_3 = (l-2*w)/((l**2-pi*a*b)*Km)
                R4 = 2*R4_1+2*R4_2+R4_3
                
                #Calculate RI
                Re = float(l /((l ** 2) * Ka))
                RI = float((1/(1 / Re - 1 / R4)-2 * R1-R2)/2)

                #Calculate Ri
                Ri = float(RI*(pi * r ** 2 / 2))
                return {"Ri": Ri}

            else:
                # Calculate thermal resistance R1.
                ex_expr=(-3*m+sqrt(-3*m**2+4*r**2))/(2*sqrt(2))
                ey_expr=(m+sqrt(-3*m**2+4*r**2))/(2*sqrt(2))
                mx_expr=(3*l-3*m-sqrt(-3*(l-m)**2+4*r**2))/(2*sqrt(2))
                my_expr=(l-m+sqrt(-3*(l-m)**2+4*r**2))/(2*sqrt(2))
                ox_expr=(3*l-3*m+sqrt(-3*(l-m)**2+4*r**2))/(2*sqrt(2))
                oy_expr=(l-m-sqrt(-3*(l-m)**2+4*r**2))/(2*sqrt(2))
                def A1_R1_1(m_val):
                    m_val = float(m_val)
                    ex_value = float(ex_expr.subs(m, m_val).evalf())
                    ey_value = float(ey_expr.subs(m, m_val).evalf())
                    def integrand(x):
                        return b / a * (a ** 2 - x ** 2) ** 0.5
                    integral_value1, _ = quad(integrand, ex_value, a)
                    return ey_value**2+2*integral_value1
                def A1_R1_2(m_val):
                    m_val = float(m_val)
                    ex_value = float(ex_expr.subs(m, m_val).evalf())
                    ey_value = float(ey_expr.subs(m, m_val).evalf())
                    ox_value = float(ox_expr.subs(m, m_val).evalf())
                    oy_value = float(oy_expr.subs(m, m_val).evalf())
                    mx_value = float(mx_expr.subs(m, m_val).evalf())
                    my_value = float(my_expr.subs(m, m_val).evalf())
                    def integrand(x):
                        return b / a * (a ** 2 - x ** 2) ** 0.5
                    integral_value2, _ = quad(integrand, ex_value, mx_value)
                    integral_value3, _ = quad(integrand, ox_value, a)
                    return ey_value**2+2*integral_value2+(oy_value+my_value)*(ox_value-mx_value)+2*integral_value3
                def integrand_R1_1(m_val):
                    try:
                        A1_val1 = float(A1_R1_1(m_val)) 
                        Vf1 = float(pi * (r**2 - m_val**2) / (4 * A1_val1))
                        Vm1 = 1 - Vf1
                        K1 = Vf1 * Kf + Vm1 * Km
                        K1 = float(K1)
                        result_1 = 1 / (A1_val1 * K1)
                        result_1 = float(result_1)
                        return result_1
                    except Exception as e:
                        print(f"Error in integrand_R1 for m_val={m_val}: {e}")
                        raise  # Re-raise the exception for debugging purposes
                def integrand_R1_2(m_val):
                    try:
                        A1_val2 = float(A1_R1_2(m_val)) 
                        Vf1 = float(pi * (r**2 - m_val**2) / (4 * A1_val2))
                        Vm1 = 1 - Vf1
                        K1 = float(Vf1 * Kf + Vm1 * Km)
                        result_2 = float(1 / (A1_val2 * K1))
                        return result_2
                    except Exception as e:
                        print(f"Error in integrand_R1 for m_val={m_val}: {e}")
                        raise  # Re-raise the exception for debugging purposes
                R1_1, _ = quad(integrand_R1_1, 0, l-w)
                R1_2, _ = quad(integrand_R1_2, l-w, l-r)
                R1 = R1_1+R1_2

                # Calculate thermal resistance R2.
                def A2_R2(m_val):
                    ex_value = float(ex_expr.subs(m, m_val).evalf())
                    mx_value = float(mx_expr.subs(m, m_val).evalf())
                    ey_value = float(ey_expr.subs(m, m_val).evalf())
                    my_value = float(my_expr.subs(m, m_val).evalf())
                    def integrand(x):
                        return b / a * (a ** 2 - x ** 2) ** 0.5
                    integral_value, _ = quad(integrand, ex_value, mx_value)
                    return ey_value**2+2 * integral_value+my_value**2
                def integrand_R2(m_val):
                    try:
                        A2_val = float(A2_R2(m_val))
                        Vf2 = float(pi * (r**2 - m_val**2+r**2-(l-m_val)**2) / (4 * A2_val))
                        Vm2 = 1 - Vf2
                        K2 = float(Vf2 * Kf + Vm2 * Km)
                        result_2 = float(1 / (A2_val * K2))
                        return result_2
                    except Exception as e:
                        raise
                R2, _ = quad(integrand_R2, l-r, r)
                # Calculate thermal resistance R4.
                def integrand_R4_1(m_val):
                    try:
                        A4_val1 = float(l**2-A1_R1_1(m_val))
                        result_R4_1 = float(1 / (A4_val1 * Km))
                        return result_R4_1
                    except Exception as e:
                        print(f"Error in integrand_R4_1 for m_val={m_val}: {e}")
                        raise
                R4_1, _ = quad(integrand_R4_1, 0, l-w)
                def integrand_R4_2(m_val):
                    try:
                        A4_val2 = float(l**2-A1_R1_2(m_val))
                        result_R4_2 = float(1 / (A4_val2 * Km))
                        return result_R4_2
                    except Exception as e:
                        print(f"Error in integrand_R4_1 for m_val={m_val}: {e}")
                        raise
                R4_2, _ = quad(integrand_R4_2, l-w, l-r)
                def integrand_R4_3(m_val):
                    try:
                        A4_val3 = float(l**2-A2_R2(m_val))
                        result_R4_3 = float(1 / (A4_val3 * Km))
                        return result_R4_3
                    except Exception as e:
                        print(f"Error in integrand_R4_1 for m_val={m_val}: {e}")
                        raise
                R4_3,_ = quad(integrand_R4_2, l-r, r)
                R4 = 2*R4_1+2*R4_2+R4_3
                #Calculate RI
                Re = float(l /((l ** 2) * Ka))
                RI = float((1/(1 / Re - 1 / R4)-2 * R1-R2)/2)

                #Calculate Ri
                Ri = float(RI*(pi * r ** 2 / 2))
                return {"Ri": Ri}
            
        return jsonify({"Ri": Ri})
    except Exception as e:
        return jsonify({"error": str(e)}), 400
        


if __name__ == "__main__":   
    # Start the Flask application.
    app.run(debug=True)

