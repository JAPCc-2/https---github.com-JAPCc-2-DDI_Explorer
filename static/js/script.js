// Menu toggle for mobile
let menu = document.querySelector('#menu-btn');
let navbar = document.querySelector('.navbar');

menu.onclick = () => {
    menu.classList.toggle('fa-times');
    navbar.classList.toggle('active');
};

// Close menu when scrolling
window.onscroll = () => {
    menu.classList.remove('fa-times');
    navbar.classList.remove('active');
};

// Interaction Checker (REAL backend connection)
function checkInteraction() {
    const drug1 = document.getElementById('drug1').value.trim();
    const drug2 = document.getElementById('drug2').value.trim();
    const result = document.getElementById('interaction-result');
    const spinner = document.getElementById('spinner');

    if (!drug1 || !drug2) {
        result.innerHTML = '<span style="color:red;">Please enter both drug names.</span>';
        spinner.style.display = 'none';
        return;
    }

    // Show spinner
    spinner.style.display = 'block';
    result.innerHTML = '';

    // Send request to Flask backend
    fetch('http://127.0.0.1:5000/predict', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json'
        },
        body: JSON.stringify({ drug1: drug1, drug2: drug2 })
    })
    .then(response => response.json())
    .then(data => {
        spinner.style.display = 'none';
        if (data.error) {
            result.innerHTML = `<span style="color:red;">${data.error}</span>`;
        } else {
            result.innerHTML = `
                <b>Result:</b> ${data.interaction}<br>
                <b>Confidence:</b> ${(data.confidence * 100).toFixed(2)}%
            `;
        }
    })
    .catch(error => {
        spinner.style.display = 'none';
        console.error('Error:', error);
        result.innerHTML = '<span style="color:red;">Error checking interaction. Please try again.</span>';
    });
}
