# Bioinformatics Matrix Prediction Application

A Python application that performs discrete step predictions using a square transition matrix and an initial state vector. This type of model is commonly used in bioinformatics for:

- Markov chain state transitions
- Population dynamics modeling
- Protein conformation predictions
- Gene expression modeling
- Epidemiological modeling

## Features

- Accepts arbitrary-sized square matrices (2×2, 3×3, 4×4, etc.)
- Performs exactly 5 discrete prediction steps
- Supports JSON input/output
- Automatically saves results to `results.json` in the ex1 folder
- Validates input dimensions and provides clear error messages
- Displays formatted results on screen
- Works from any directory (but recommended to run from ex1 folder)

## Requirements

- Python 3.6+
- NumPy

Install dependencies:
```bash
pip install numpy
```

## Usage

**IMPORTANT:** Always run the script from within the `ex1` folder:
```bash
cd ex1
python predictor.py [options]
```

### Basic Usage (Default Example)

Run with built-in example data:
```bash
python predictor.py
```

This will automatically create `results.json` in the ex1 folder with the prediction results.

### Using Custom Input File

Provide a JSON file with your matrix and initial vector:
```bash
python predictor.py example_input.json
```

### Custom Output File (Optional)

By default, results are saved to `results.json` in the ex1 folder. You can specify a custom output file:
```bash
python predictor.py example_input.json my_output.json
```

## Input Format

Create a JSON file with the following structure:

```json
{
  "matrix": [
    [0.8, 0.15, 0.05],
    [0.1, 0.7, 0.2],
    [0.1, 0.15, 0.75]
  ],
  "initial_vector": [100, 0, 0]
}
```

- **matrix**: Square transition matrix (n × n)
- **initial_vector**: Initial state vector (n elements)

## Example Files

The ex1 folder includes several example input files:

### example_input.json
Demonstrates a 3×3 transition matrix with initial population of 100 in state 0.

### example_markov.json
Demonstrates a Markov chain with probability distribution starting at state 1.

### test_2x2.json
Simple 2×2 matrix example for testing.

### test_4x4.json
4×4 matrix example demonstrating arbitrary size support.

## Mathematical Model

At each step, the new state is calculated as:

```
State(t+1) = Matrix × State(t)
```

The application performs this calculation 5 times, showing the state evolution from step 0 (initial) to step 5.

## Examples

All commands assume you are in the `ex1` folder:

### Running with example_input.json:
```bash
python predictor.py example_input.json
```

### Running with example_markov.json:
```bash
python predictor.py example_markov.json
```

### Running with 2×2 test case:
```bash
python predictor.py test_2x2.json
```

### Running with 4×4 test case:
```bash
python predictor.py test_4x4.json
```

## Output

The application automatically:
1. Displays the transition matrix on screen
2. Shows initial state (Step 0)
3. Shows predicted states for Steps 1-5
4. **Saves complete results to `results.json`** in the ex1 folder

### results.json Format

The output JSON file contains:
```json
{
  "matrix": [...],
  "initial_vector": [...],
  "predictions": [
    [...],  // Step 0 (initial)
    [...],  // Step 1
    [...],  // Step 2
    [...],  // Step 3
    [...],  // Step 4
    [...]   // Step 5
  ]
}
```

## Bioinformatics Applications

This matrix prediction model can be applied to:

1. **Protein Folding**: Model transitions between conformational states
2. **Population Genetics**: Predict allele frequency changes over generations
3. **Epidemic Modeling**: Track disease state distributions (susceptible/infected/recovered)
4. **Gene Regulatory Networks**: Model gene expression states
5. **Sequence Evolution**: Predict nucleotide/amino acid substitutions

## Implementation Details

- Uses NumPy for efficient matrix operations
- Validates matrix is square and dimensions match
- Handles both row and column vectors
- Supports floating-point precision
- Error handling for invalid inputs
- Automatic path resolution for input/output files
- All results automatically saved to ex1 folder

## Quick Start

1. Navigate to the ex1 folder:
   ```bash
   cd C:\Users\alina\Desktop\L12\ex1
   ```

2. Run with an example file:
   ```bash
   python predictor.py example_input.json
   ```

3. Check the results in `results.json`

## Troubleshooting

**Error: "can't open file"**
- Make sure you are in the ex1 folder before running the script
- Use `cd ex1` to navigate to the correct folder
- Don't include `ex1/` in the command when you're already in the ex1 folder

**Error: "FileNotFoundError"**
- Verify the input JSON file exists in the ex1 folder
- Check the filename spelling and extension (.json)
