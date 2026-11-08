"""
server.py
-------------------------------------------------------------------
Production REST API for Noncoding Cancer-Driver Mutation Hotspots.
Serves ML predictions (Melanoma vs. Breast Cancer).
-------------------------------------------------------------------
"""

import os
from pathlib import Path
import pandas as pd
import numpy as np
from flask import Flask, request, jsonify
from flask_cors import CORS
from sklearn.ensemble import RandomForestClassifier
import joblib

app = Flask(__name__)
# Enable CORS for external frontend calls
CORS(app)

BASE_DIR = Path(__file__).resolve().parent

# Flexible path lookup for matrix CSV
POSSIBLE_MATRIX_PATHS = [
    BASE_DIR / "ml_hotspot_training_matrix.csv",
    BASE_DIR / "results" / "tables" / "ml_hotspot_training_matrix.csv",
    Path("ml_hotspot_training_matrix.csv")
]

MODEL_FILE = BASE_DIR / "rf_hotspot_classifier.pkl"

FEATURE_COLS = [
    "Mutated_Donors", 
    "Fold_Enrichment", 
    "Absolute_Distance_to_TSS", 
    "Is_Promoter", 
    "Is_Enhancer", 
    "Has_TF_Impact", 
    "Motif_Gain", 
    "Motif_Loss"
]

MODEL = None


def find_matrix_file():
    for p in POSSIBLE_MATRIX_PATHS:
        if p.exists():
            return p
    return None


def initialize_engine():
    """Loads pre-trained model or trains a fresh Random Forest classifier on boot."""
    global MODEL
    matrix_path = find_matrix_file()

    if MODEL_FILE.exists():
        print(f"[*] Loading serialized model from disk: {MODEL_FILE}")
        MODEL = joblib.load(MODEL_FILE)
    elif matrix_path is not None:
        print(f"[*] Training Random Forest Classifier from: {matrix_path}")
        df = pd.read_csv(matrix_path)
        X = df[FEATURE_COLS]
        y = df["Cancer_Type"]

        clf = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42)
        clf.fit(X, y)
        
        joblib.dump(clf, MODEL_FILE)
        MODEL = clf
        print("[+] Model successfully trained and serialized.")
    else:
        print("[!] Warning: Matrix CSV not found. API initialized in standby mode.")


# Initialize model on startup
initialize_engine()


@app.route('/', methods=['GET'])
def health_check():
    return jsonify({
        "status": "online",
        "service": "Noncoding Driver Mutation Classification API",
        "model_loaded": MODEL is not None
    }), 200


@app.route('/predict', methods=['POST'])
def predict():
    if MODEL is None:
        return jsonify({
            "status": "error",
            "message": "Engine unavailable. Training matrix CSV was not found on server."
        }), 500

    try:
        data = request.get_json(force=True)
        if not data:
            return jsonify({"status": "error", "message": "Empty JSON payload."}), 400

        missing_fields = [col for col in FEATURE_COLS if col not in data]
        if missing_fields:
            return jsonify({
                "status": "error",
                "message": f"Missing required features: {missing_fields}"
            }), 400

        input_data = {col: [float(data[col])] for col in FEATURE_COLS}
        input_df = pd.DataFrame(input_data)

        prediction = MODEL.predict(input_df)[0]
        probabilities = MODEL.predict_proba(input_df)[0]
        class_names = MODEL.classes_

        prob_distribution = {
            cls: float(np.round(prob * 100, 2)) 
            for cls, prob in zip(class_names, probabilities)
        }
        confidence = float(np.round(max(probabilities) * 100, 2))

        return jsonify({
            "status": "success",
            "predicted_cancer_type": prediction,
            "confidence_score": confidence,
            "probability_breakdown": prob_distribution
        }), 200

    except Exception as e:
        return jsonify({"status": "error", "message": str(e)}), 400


if __name__ == '__main__':
    port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=port, debug=False)