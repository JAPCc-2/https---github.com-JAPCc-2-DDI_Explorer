Introduction

DDI Explorer helps medical professionals and researchers quickly assess potential drug-drug interactions based on molecular structure data. The backend of the system is powered by Python, Flask (for the web framework), Random Forest (for the machine learning model), and RDKit (for chemical data processing).

Prerequisites

To get started, make sure your system meets the following requirements:

Software Requirements:

Python 3.8+

pip (Python Package Installer)
To install required Python packages.

Git

To clone the repository from GitHub. Install from here.

Virtual Environment (Optional but recommended)

It’s highly recommended to set up a virtual environment to avoid conflicts between dependencies.



Setting Up the Project


Step 1: Clone the Repository
First, clone the project repository to your local machine using the following command:


bash
Copy
git clone https---github.com-JAPCc-2-DDI_Explorer.git



Navigate into the project folder:


bash
Copy
cd DDI_Explorer


Step 2: Set Up a Virtual Environment (Optional)
Install virtualenv if it’s not already installed:

bash
Copy
pip install virtualenv


Create a new virtual environment:

bash
Copy
python3 -m venv venv


Activate the virtual environment:

For macOS/Linux:

bash
Copy
source venv/bin/activate


For Windows:

bash
Copy
venv\Scripts\activate


Step 3: Install Required Dependencies
Install the necessary Python packages using pip:

bash
Copy
pip install -r requirements.txt


This will install all dependencies listed in the requirements.txt file (Flask, scikit-learn, RDKit, requests, etc.).

Running the Application

Step 1: Start the Flask Development Server
Inside your project directory, run the following command to start the Flask server:

bash
Copy
python app.py


You should see an output similar to this:

csharp
Copy
* Running on http://127.0.0.1:5000 (Press CTRL+C to quit)
* Restarting with stat
This means the server is running locally on localhost port 5000.

Open your browser and go to
http://127.0.0.1:5000
You should see the DDI Explorer web app!

Step 2: Interacting with the Web App
Once the web app is open in your browser, enter the names of two drugs in the input fields, such as:

Drug 1: Atorvastatin

Drug 2: Simvastatin

Click the "Analyze Interaction" button.

After a brief delay, the app will display whether there is a potential interaction between the two drugs, along with a confidence score.
