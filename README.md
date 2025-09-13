# Database Normalization Tool

A comprehensive desktop application for database normalization analysis and schema design. This tool helps database designers and students understand and apply normalization principles from First Normal Form (1NF) through Boyce-Codd Normal Form (BCNF).

## Features

- **Attribute Closure Computation**: Calculate the closure of individual attributes and attribute sets
- **Candidate Key Identification**: Automatically find all candidate keys for a given schema
- **Complete Normalization Process**: Step-by-step normalization from 1NF to BCNF
- **Interactive GUI**: User-friendly interface built with Tkinter
- **Detailed Analysis**: Comprehensive explanations of each normalization step
  
## Snapshots
![image_alt](https://github.com/hardattmangrola/Database-Normalization-Tool/blob/be3db349bd88b019b91919dbd60e83dfc90293aa/ui1.png)
![image_alt](https://github.com/hardattmangrola/Database-Normalization-Tool/blob/be3db349bd88b019b91919dbd60e83dfc90293aa/ui2.png)

## Requirements

- Python 3.6 or higher
- Tkinter (included with most Python installations)

## Installation

1. Clone or download this repository
2. Ensure Python 3.6+ is installed on your system
3. No additional dependencies required

## Usage

### Running the Application

```bash
python main.py
```

### Input Format

1. **Schema Definition**: Enter comma-separated attributes (e.g., `A, B, C, D, E`)

2. **Functional Dependencies**: Enter one dependency per line using the format:
   - Single attribute: `A -> B`
   - Multiple attributes: `A,B -> C,D`

### Example Input

**Schema:**
```
A, B, C, D, E
```

**Functional Dependencies:**
```
A -> B
B -> C
D -> E
A,D -> C
```

### Available Operations

- **Find Attribute Closures**: Computes and displays the closure of each attribute and attribute combination
- **Find Candidate Keys**: Identifies all candidate keys using systematic analysis
- **Normalize Schema**: Performs complete normalization process with detailed explanations
- **Clear Results**: Clears the output area for new analysis

## Normalization Process

The tool performs normalization through the following steps:

1. **First Normal Form (1NF)**: Ensures atomic attributes
2. **Second Normal Form (2NF)**: Eliminates partial dependencies
3. **Third Normal Form (3NF)**: Removes transitive dependencies
4. **Boyce-Codd Normal Form (BCNF)**: Ensures every determinant is a superkey

Each step includes:
- Detailed explanations of the normal form requirements
- Identification of violations
- Step-by-step decomposition process
- Final normalized schema

## Technical Implementation

- **GUI Framework**: Tkinter for cross-platform compatibility
- **Algorithm**: Implements standard database normalization algorithms
- **Data Structures**: Uses sets and frozensets for efficient attribute manipulation
- **Error Handling**: Comprehensive input validation and error reporting

## Educational Value

This tool is particularly useful for:
- Database design students learning normalization concepts
- Database administrators reviewing schema designs
- Software developers understanding database optimization
- Academic instructors teaching database theory

## License

This project is open source and available under standard licensing terms.
