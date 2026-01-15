"""
Bioinformatics Matrix Prediction Application
Performs discrete step predictions using a transition matrix and initial state vector.
"""

import numpy as np
import json
import sys
import os


class MatrixPredictor:
    """
    A predictor that applies matrix multiplication for discrete step predictions.
    Common in bioinformatics for Markov chain models, state transitions, etc.
    """

    def __init__(self, matrix, initial_vector):
        """
        Initialize the predictor with a transition matrix and initial state vector.

        Args:
            matrix: Square transition matrix (n x n)
            initial_vector: Initial state vector (n x 1)
        """
        self.matrix = np.array(matrix, dtype=float)
        self.initial_vector = np.array(initial_vector, dtype=float)

        # Validate inputs
        self._validate_inputs()

    def _validate_inputs(self):
        """Validate that matrix is square and dimensions match."""
        if self.matrix.ndim != 2:
            raise ValueError("Matrix must be 2-dimensional")

        if self.matrix.shape[0] != self.matrix.shape[1]:
            raise ValueError(f"Matrix must be square, got shape {self.matrix.shape}")

        if self.initial_vector.ndim == 1:
            self.initial_vector = self.initial_vector.reshape(-1, 1)

        if self.matrix.shape[1] != self.initial_vector.shape[0]:
            raise ValueError(
                f"Matrix columns ({self.matrix.shape[1]}) must match "
                f"vector rows ({self.initial_vector.shape[0]})"
            )

    def predict_step(self, current_vector):
        """
        Perform a single prediction step: new_state = Matrix × current_state

        Args:
            current_vector: Current state vector

        Returns:
            Next state vector after applying the transition matrix
        """
        return np.dot(self.matrix, current_vector)

    def predict_n_steps(self, n_steps=5):
        """
        Perform n discrete prediction steps.

        Args:
            n_steps: Number of steps to predict (default: 5)

        Returns:
            List of state vectors at each step (including initial state)
        """
        states = [self.initial_vector.copy()]
        current_state = self.initial_vector.copy()

        for step in range(n_steps):
            current_state = self.predict_step(current_state)
            states.append(current_state.copy())

        return states

    def display_predictions(self, states):
        """
        Display the prediction results in a formatted manner.

        Args:
            states: List of state vectors at each step
        """
        print("=" * 60)
        print("BIOINFORMATICS MATRIX PREDICTION RESULTS")
        print("=" * 60)
        print(f"\nTransition Matrix ({self.matrix.shape[0]}x{self.matrix.shape[1]}):")
        print(self.matrix)
        print("\n" + "=" * 60)

        for i, state in enumerate(states):
            if i == 0:
                print(f"\nInitial State (Step 0):")
            else:
                print(f"\nPrediction Step {i}:")
            print(state.flatten())
            print("-" * 60)

        print("\n" + "=" * 60)
        print("PREDICTION COMPLETE")
        print("=" * 60)


def load_data_from_json(filename):
    """
    Load matrix and initial vector from JSON file.

    Expected JSON format:
    {
        "matrix": [[...], [...], ...],
        "initial_vector": [...]
    }
    """
    with open(filename, 'r') as f:
        data = json.load(f)

    return data['matrix'], data['initial_vector']


def main():
    """Main application entry point."""
    # Get the directory where this script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))

    if len(sys.argv) > 1:
        # Load from JSON file if provided
        input_file = sys.argv[1]
        # Try the file as-is first, then relative to script directory
        if not os.path.isabs(input_file):
            if os.path.exists(input_file):
                input_file = os.path.abspath(input_file)
            else:
                input_file = os.path.join(script_dir, input_file)
        print(f"Loading data from {input_file}...")
        matrix, initial_vector = load_data_from_json(input_file)
    else:
        # Use default example data
        print("No input file provided. Using default example data.")
        print("Usage: python predictor.py [input_file.json]\n")

        # Example: Simple 3-state transition matrix (e.g., protein conformations)
        matrix = [
            [0.7, 0.2, 0.1],
            [0.1, 0.6, 0.3],
            [0.2, 0.2, 0.6]
        ]

        # Initial state vector (e.g., initial population distribution)
        initial_vector = [1.0, 0.0, 0.0]

    # Create predictor and run predictions
    predictor = MatrixPredictor(matrix, initial_vector)
    states = predictor.predict_n_steps(n_steps=5)
    predictor.display_predictions(states)

    # Optional: Save results
    if len(sys.argv) > 2:
        output_file = sys.argv[2]
        # Try the file as-is first, then relative to script directory
        if not os.path.isabs(output_file):
            output_dir = os.path.dirname(output_file)
            if output_dir and os.path.exists(output_dir):
                output_file = os.path.abspath(output_file)
            else:
                output_file = os.path.join(script_dir, output_file)
    else:
        # Default output file in script directory
        output_file = os.path.join(script_dir, 'results.json')

    # Always save results to results.json
    results = {
        'matrix': matrix if isinstance(matrix, list) else matrix.tolist(),
        'initial_vector': initial_vector if isinstance(initial_vector, list) else initial_vector.tolist(),
        'predictions': [state.flatten().tolist() for state in states]
    }
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {output_file}")


if __name__ == "__main__":
    main()
