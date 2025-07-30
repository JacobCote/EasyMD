# Testing Framework

## Overview
EasyMD includes a comprehensive testing framework that ensures reliability, correctness, and robustness across all components. The test suite covers unit tests, integration tests, validation tests, and edge case handling, providing confidence in the software's behavior under various conditions.

## Testing Architecture

### Test Organization
```
tests/
├── data/                              # Test data files
├── test_analysis.py                   # Analysis component tests
├── test_arg_manager.py               # Argument manager tests
├── test_data_utils.py                # Data utility tests
├── test_info.py                      # Info component tests
├── test_integration.py               # Integration tests
├── test_missing_residue_handling.py  # Missing residue tests
├── test_runners.py                   # Runner component tests
├── test_sim_runner.py                # SimRunner tests
├── test_sys_generator.py             # SysGenerator tests
├── test_terminal_residue_limit.py    # Terminal residue tests
├── test_utils.py                     # Utility function tests
├── test_validation_display.py        # Validation display tests
├── test_validation_warnings_errors.py # Error handling tests
├── test_validation.py               # General validation tests
├── test_water_handling.py           # Water handling tests
├── test_your_structure.py           # Structure validation tests
├── testArgManager.py                # Legacy ArgManager tests
├── testSimRunner.py                 # Legacy SimRunner tests
├── testSysGenerator.py              # Legacy SysGenerator tests
└── debug_terminal_residues.py       # Debug utilities
```

### Test Categories

#### 1. Unit Tests
- **Component Isolation**: Test individual components in isolation
- **Function-Level Testing**: Validate specific function behavior
- **Mock Dependencies**: Use mock objects to isolate functionality
- **Edge Case Coverage**: Test boundary conditions and error cases

#### 2. Integration Tests
- **Component Interaction**: Test interaction between components
- **Workflow Testing**: Validate complete simulation workflows
- **End-to-End Testing**: Test from input to output
- **Real Data Testing**: Use actual PDB files and configurations

#### 3. Validation Tests
- **Input Validation**: Test argument parsing and validation
- **Error Handling**: Verify proper error detection and reporting
- **Warning Systems**: Test warning generation and display
- **Configuration Testing**: Validate configuration file handling

#### 4. Regression Tests
- **Behavior Preservation**: Ensure changes don't break existing functionality
- **Output Consistency**: Verify consistent output across versions
- **Performance Testing**: Monitor performance characteristics
- **Compatibility Testing**: Test with different dependencies

## Component-Specific Testing

### 1. ArgManager Testing (test_arg_manager.py)
```python
class TestArgManager:
    """Test suite for argument parsing and validation."""
    
    def test_valid_configuration(self):
        """Test that valid configurations are accepted."""
        # Test various valid parameter combinations
        
    def test_invalid_configuration(self):
        """Test that invalid configurations are rejected."""
        # Test missing required parameters
        # Test conflicting parameters
        # Test out-of-range values
        
    def test_config_file_loading(self):
        """Test YAML configuration file loading."""
        # Test valid config files
        # Test invalid YAML syntax
        # Test missing config files
        
    def test_parameter_validation(self):
        """Test individual parameter validation."""
        # Test temperature ranges
        # Test step size validation
        # Test file path validation
```

**Test Coverage:**
- **Parameter Validation**: All argument types and ranges
- **Configuration Files**: YAML loading and parsing
- **Error Messages**: Clear, actionable error reporting
- **Edge Cases**: Boundary conditions and invalid inputs

### 2. SysGenerator Testing (test_sys_generator.py)
```python
class TestSysGenerator:
    """Test suite for system preparation and generation."""
    
    def test_protein_preparation(self):
        """Test protein-only system preparation."""
        # Test structure fixing
        # Test missing residue handling
        # Test solvation
        
    def test_complex_preparation(self):
        """Test protein-ligand complex preparation."""
        # Test ligand extraction
        # Test charge calculation
        # Test system combination
        
    def test_missing_residue_handling(self):
        """Test missing residue detection and handling."""
        # Test PDBFixer integration
        # Test different strategies
        # Test terminal residue limits
        
    def test_coordinate_validation(self):
        """Test coordinate validation and error detection."""
        # Test NaN detection
        # Test infinite coordinates
        # Test reasonable ranges
```

**Test Coverage:**
- **Structure Preparation**: Protein and complex systems
- **Missing Residues**: PDBFixer integration and strategies
- **Coordinate Validation**: NaN and infinite value detection
- **Force Field Assignment**: Proper parameterization

### 3. Info Component Testing (test_info.py)
```python
class TestInfoRunner:
    """Test suite for PDB structure analysis."""
    
    def test_pdb_parsing(self):
        """Test PDB file parsing and data extraction."""
        # Test standard PDB files
        # Test compressed files
        # Test malformed files
        
    def test_missing_residue_detection(self):
        """Test missing residue detection accuracy."""
        # Test PDBFixer integration
        # Test gap detection fallback
        # Test chain-specific analysis
        
    def test_ligand_identification(self):
        """Test ligand and small molecule identification."""
        # Test ligand detection
        # Test solvent filtering
        # Test multi-chain ligands
        
    def test_output_formatting(self):
        """Test output generation and formatting."""
        # Test detailed reports
        # Test summary format
        # Test JSON output
```

**Test Coverage:**
- **PDB Parsing**: Various file formats and structures
- **Analysis Accuracy**: Comparison with known results
- **Output Formats**: All supported output types
- **Error Handling**: Malformed or incomplete files

### 4. Analysis Testing (test_analysis.py)
```python
class TestAnalysisRunner:
    """Test suite for trajectory analysis."""
    
    def test_trajectory_loading(self):
        """Test trajectory file loading and combination."""
        # Test single trajectory files
        # Test multiple trajectory combination
        # Test topology loading
        
    def test_rmsd_calculation(self):
        """Test RMSD calculation accuracy."""
        # Test different atom selections
        # Test reference frame selection
        # Test known RMSD values
        
    def test_rmsf_calculation(self):
        """Test RMSF calculation accuracy."""
        # Test per-residue calculations
        # Test CA atom selection
        # Test statistical measures
        
    def test_output_generation(self):
        """Test plot and data file generation."""
        # Test plot creation
        # Test data export
        # Test format options
```

**Test Coverage:**
- **Trajectory Processing**: Loading and manipulation
- **Analysis Accuracy**: Numerical correctness
- **Output Generation**: Plots and data files
- **Performance**: Large trajectory handling

### 5. Validation Testing (test_validation.py)
```python
class TestValidation:
    """Test suite for input validation systems."""
    
    def test_error_detection(self):
        """Test error detection and reporting."""
        # Test various error conditions
        # Test error message clarity
        # Test error categorization
        
    def test_warning_system(self):
        """Test warning generation and display."""
        # Test warning conditions
        # Test warning vs error distinction
        # Test warning suppression
        
    def test_validation_display(self):
        """Test validation result display."""
        # Test success messages
        # Test error formatting
        # Test color output
```

**Test Coverage:**
- **Error Detection**: Comprehensive error identification
- **Warning Systems**: Appropriate warning generation
- **User Experience**: Clear, helpful messages
- **Display Formatting**: Proper output formatting

## Testing Utilities and Infrastructure

### 1. Test Data Management
```python
# Test data organization
tests/data/
├── test_protein.pdb          # Simple protein structure
├── test_complex.pdb          # Protein-ligand complex
├── malformed.pdb            # Intentionally malformed file
├── missing_residues.pdb     # Structure with missing residues
├── config_valid.yml         # Valid configuration file
├── config_invalid.yml       # Invalid configuration file
└── trajectory_files/        # Test trajectory data
```

### 2. Mock Objects and Fixtures
```python
@pytest.fixture
def mock_config():
    """Create mock configuration for testing."""
    config = Mock()
    config.protein = 'tests/data/test_protein.pdb'
    config.steps = 1000
    config.temperature = 300
    config.solvate = True
    return config

@pytest.fixture
def test_pdb_file():
    """Provide test PDB file path."""
    return 'tests/data/test_protein.pdb'
```

### 3. Test Utilities
```python
def assert_file_exists(filepath):
    """Assert that a file exists and is readable."""
    assert os.path.exists(filepath), f"File does not exist: {filepath}"
    assert os.path.isfile(filepath), f"Path is not a file: {filepath}"

def assert_valid_pdb(filepath):
    """Assert that a file is a valid PDB file."""
    with open(filepath, 'r') as f:
        content = f.read()
        assert 'ATOM' in content or 'HETATM' in content
        assert content.strip().endswith('END')

def compare_numerical_results(result1, result2, tolerance=1e-6):
    """Compare numerical results with tolerance."""
    assert abs(result1 - result2) < tolerance
```

## Test Execution and Coverage

### 1. Running Tests
```bash
# Run all tests
pytest src/EasyMD/tests/

# Run specific test file
pytest src/EasyMD/tests/test_arg_manager.py

# Run with coverage
pytest --cov=src/EasyMD src/EasyMD/tests/

# Run with verbose output
pytest -v src/EasyMD/tests/

# Run specific test function
pytest src/EasyMD/tests/test_info.py::TestInfoRunner::test_missing_residue_detection
```

### 2. Coverage Analysis
```bash
# Generate coverage report
pytest --cov=src/EasyMD --cov-report=html src/EasyMD/tests/

# Coverage targets
# - Overall coverage: >90%
# - Critical components: >95%
# - Error handling: 100%
```

### 3. Continuous Integration
```yaml
# Example CI configuration
name: Test Suite
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - name: Set up Python
        uses: actions/setup-python@v2
        with:
          python-version: 3.8
      - name: Install dependencies
        run: pip install -r requirements.txt
      - name: Run tests
        run: pytest src/EasyMD/tests/
```

## Quality Assurance Measures

### 1. Test Quality Standards
- **Comprehensive Coverage**: All public functions tested
- **Edge Case Testing**: Boundary conditions and error cases
- **Integration Testing**: Component interaction validation
- **Performance Testing**: Resource usage and timing

### 2. Test Maintenance
- **Regular Updates**: Tests updated with code changes
- **Regression Prevention**: New tests for bug fixes
- **Documentation**: Clear test descriptions and purposes
- **Cleanup**: Removal of obsolete or redundant tests

### 3. Error Simulation
```python
def test_nan_coordinate_handling():
    """Test handling of NaN coordinates."""
    # Create structure with NaN coordinates
    # Verify proper error detection
    # Check error message quality
    
def test_missing_file_handling():
    """Test handling of missing input files."""
    # Test with non-existent files
    # Verify appropriate error messages
    # Check graceful failure
```

## Testing Best Practices

### 1. Test Design Principles
- **Isolation**: Tests don't depend on each other
- **Repeatability**: Tests produce consistent results
- **Clarity**: Test purpose is immediately clear
- **Maintainability**: Tests are easy to update and modify

### 2. Mock Usage
```python
@patch('EasyMD.utils.get_platform')
def test_platform_selection(mock_get_platform):
    """Test platform selection with mocked platform."""
    mock_platform = Mock()
    mock_platform.getName.return_value = 'CUDA'
    mock_get_platform.return_value = mock_platform
    
    # Test platform-dependent functionality
```

### 3. Parameterized Testing
```python
@pytest.mark.parametrize("temperature,expected", [
    (300, True),   # Valid temperature
    (0, False),    # Invalid temperature
    (-10, False),  # Negative temperature
    (1000, False), # Extremely high temperature
])
def test_temperature_validation(temperature, expected):
    """Test temperature validation with various values."""
    result = validate_temperature(temperature)
    assert result == expected
```

## Benefits of Comprehensive Testing

1. **Reliability**: Ensures software behaves correctly under various conditions
2. **Regression Prevention**: Catches bugs introduced by changes
3. **Documentation**: Tests serve as executable documentation
4. **Confidence**: Provides confidence in software quality
5. **Maintainability**: Makes refactoring safer and easier
6. **User Experience**: Ensures consistent, predictable behavior
7. **Development Speed**: Catches issues early in development

## Test Results and Validation

### Current Test Status
- **Total Tests**: 21 test files covering all major components
- **Coverage**: Comprehensive coverage of core functionality
- **Integration**: Full workflow testing from input to output
- **Validation**: Extensive input validation and error handling

### Key Achievements
- **PDBFixer Integration**: Validated accurate missing residue detection
- **Error Handling**: Comprehensive error detection and user guidance
- **Component Integration**: Verified proper interaction between components
- **Real-World Testing**: Validated with actual PDB structures (4zgm.pdb, etc.)

The testing framework ensures that EasyMD maintains high quality, reliability, and user-friendliness across all its components, providing researchers with confidence in their molecular dynamics simulation results.