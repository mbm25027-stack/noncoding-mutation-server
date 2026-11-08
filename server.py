"""
server.py
-------------------------------------------------------------------
Production REST API for Noncoding Cancer-Driver Mutation Hotspots.
Serves Machine Learning predictions for external web frontends.
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
# Enable CORS so your website (e.g., https://mutationpredict.org) can talk to this API
CORS(app)

# Path setup
BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "results"
MATRIX_FILE = RESULTS_DIR / "tables" / "ml_hotspot_training_matrix.csv"
MODEL_FILE = RESULTS_DIR / "models" / "rf_hotspot_classifier.pkl"

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


def initialize_engine():
    """Loads existing pre-trained model or trains a fresh instance on startup."""
    global MODEL
    MODEL_FILE.parent.mkdir(parents=True, exist_ok=True)

    if MODEL_FILE.exists():
        print(f"[*] Loading serialized model from: {MODEL_FILE.name}")
        MODEL = joblib.load(MODEL_FILE)
    elif MATRIX_FILE.exists():
        print(f"[*] Training new Random Forest Classifier from: {MATRIX_FILE.name}")
        df = pd.read_csv(MATRIX_FILE)
        X = df[FEATURE_COLS]
        y = df["Cancer_Type"]

        clf = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42)
        clf.fit(X, y)
        
        # Save model to disk
        joblib.dump(clf, MODEL_FILE)
        MODEL = clf
        print(f"[+] Model saved to: {MODEL_FILE}")
    else:
        print("[!] Warning: Neither pre-trained model nor training matrix CSV found.")


# Initialize model on boot
initialize_engine()


@app.route('/', methods=['GET'])
def health_check():
    """Health check endpoint to verify server status on cloud hosts."""
    return jsonify({
        "status": "online",
        "service": "Noncoding Mutation Driver Prediction API",
        "model_loaded": MODEL is not None,
        "supported_features": FEATURE_COLS
    }), 200


@app.route('/predict', methods=['POST'])
def predict():
    """
    Predicts cancer cohort based on noncoding mutation features.
    Accepts JSON body matching FEATURE_COLS.
    """
    if MODEL is None:
        return jsonify({
            "status": "error",
            "message": "Engine unavailable. Training matrix not loaded."
        }), 500

    try:
        data = request.get_json(force=True)
        if not data:
            return jsonify({"status": "error", "message": "Empty JSON payload."}), 400

        # Validate input features and check for missing fields
        missing_fields = [col for col in FEATURE_COLS if col not in data]
        if missing_fields:
            return jsonify({
                "status": "error",
                "message": f"Missing required features: {missing_fields}"
            }), 400

        # Build feature vector DataFrame
        input_data = {col: [float(data[col])] for col in FEATURE_COLS}
        input_df = pd.DataFrame(input_data)

        # Inference
        prediction = MODEL.predict(input_df)[0]
        probabilities = MODEL.predict_proba(input_df)[0]
        class_names = MODEL.classes_

        # Format probability distribution
        prob_distribution = {
            cls: float(np.round(prob * 100, 2)) 
            for cls, prob in zip(class_names, probabilities)
        }
        confidence = float(np.round(max(probabilities) * 100, 2))

        return jsonify({
            "status": "success",
            "predicted_cancer_type": prediction,
            "confidence_score": confidence,
            "probability_breakdown": prob_distribution,
            "analyzed_features": data
        }), 200

    except Exception as e:
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 400


if __name__ == '__main__':
    # PORT defaults to 5000 locally, but pulls environment variable for cloud deployment
    port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=port, debug=False)