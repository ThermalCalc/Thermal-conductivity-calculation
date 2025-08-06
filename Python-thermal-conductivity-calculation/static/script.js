// Single Data Set Calculation via Form Submission
document.getElementById('single-form').addEventListener('submit', async (e) => {
    e.preventDefault();

    const formData = new FormData(e.target);
    const data = {
        r: parseFloat(formData.get('r')),
        Vf: parseFloat(formData.get('Vf')),
        Kf: parseFloat(formData.get('Kf')),
        Km: parseFloat(formData.get('Km')),
        Ri: parseFloat(formData.get('Ri'))
    };

    const response = await fetch('/calculate-single', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(data)
    });

    const responseText = await response.text();
    console.log('Response Text:', responseText);

    try {
        const result = JSON.parse(responseText);
        document.getElementById('single-result').innerHTML = `
            <h3>Calculation Results</h3>
            <p>${JSON.stringify(result)}</p>
        `;
    } catch (error) {
        console.error('Error parsing JSON:', error);
    }
});

// Batch Data Calculation Form Submission
document.getElementById('multiple-form').addEventListener('submit', async (e) => {
    e.preventDefault();

    const formData = new FormData(e.target);

    const response = await fetch('/calculate-multiple', {
        method: 'POST',
        body: formData
    });

    const result = await response.json();
    document.getElementById('multiple-result').innerHTML = `
        <h3>Calculation Results</h3>
        <p>${JSON.stringify(result)}</p>
        <a href="${result.result.filePath}" download>Download result file</a>
    `;
});

// Frontend JS for interface resistance calculation form
window.onload = function() {
    document.getElementById('interface-form').addEventListener('submit', async (e) => {
        e.preventDefault();

        const formData = new FormData(e.target);
        const response = await fetch('/calculate-interface', {
            method: 'POST',
            body: new URLSearchParams(formData), // Use URLSearchParams for form data submission
        });

        if (response.ok) {
            const result = await response.json();
            document.getElementById('interface-result').innerHTML = `
                <h3>Calculation Result</h3>
                <p>${JSON.stringify(result)}</p>
            `;
        } else {
            document.getElementById('interface-result').innerHTML = `
                <h3>Error: Could not calculate interface resistance</h3>
            `;
        }
    });
};
