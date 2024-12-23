# Contributing to CurvedShellAnalysis

## Table of Contents
1. [Code of Conduct](#code-of-conduct)
2. [Getting Started](#getting-started)
3. [Development Process](#development-process)
4. [Code Style Guide](#code-style-guide)
5. [Testing Guidelines](#testing-guidelines)
6. [Documentation](#documentation)
7. [Pull Request Process](#pull-request-process)

## Code of Conduct

### Our Pledge
We are committed to providing a friendly, safe, and welcoming environment for all contributors.

### Our Standards
- Be respectful and inclusive
- Accept constructive criticism gracefully
- Focus on what's best for the community
- Show empathy towards others

## Getting Started

### Prerequisites
1. MATLAB R2020b or later
2. Git for version control
3. GitHub account

### Setting Up Development Environment
```bash
# Clone repository
git clone https://github.com/yourusername/CurvedShellAnalysis.git
cd CurvedShellAnalysis

# Create branch
git checkout -b feature/your-feature-name
```

### Running Tests
```matlab
% Run all tests
runtests('CurvedShellAnalysis');

% Run specific test suite
runtests('CurvedShellAnalysis.Tests.ModalAnalysisTest');
```

## Development Process

### Branch Naming Convention
- `feature/`: New features
- `bugfix/`: Bug fixes
- `docs/`: Documentation updates
- `test/`: Test additions or modifications
- `refactor/`: Code refactoring

### Commit Messages
Follow the conventional commits specification:
```
type(scope): description

[optional body]

[optional footer]
```

Types:
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation
- `style`: Formatting
- `refactor`: Code restructuring
- `test`: Tests
- `chore`: Maintenance

Example:
```
feat(modal): add damping to modal analysis

- Implemented Rayleigh damping
- Added damping ratio parameter
- Updated documentation

Closes #123
```

## Code Style Guide

### File Organization
```matlab
% Class definition file structure
classdef ClassName < handle
    % Class description
    %
    % Properties:
    %   prop1 - Description
    %   prop2 - Description
    %
    % Methods:
    %   method1 - Description
    %   method2 - Description
    
    properties
        % Property documentation
        Property1
        Property2
    end
    
    methods
        function obj = ClassName(param1, param2)
            % Constructor
            % 
            % Parameters:
            %   param1 - Description
            %   param2 - Description
        end
        
        function output = method1(obj, input)
            % Method description
            %
            % Parameters:
            %   input - Description
            %
            % Returns:
            %   output - Description
        end
    end
end
```

### Naming Conventions
1. Classes: PascalCase
   ```matlab
   SphericalSurface
   ModalAnalysis
   ```

2. Functions/Methods: camelCase
   ```matlab
   calculateStiffness
   plotResults
   ```

3. Variables: camelCase
   ```matlab
   nodeCoordinates
   elementMatrix
   ```

4. Constants: UPPER_CASE
   ```matlab
   MAX_ITERATIONS
   TOLERANCE
   ```

### Code Formatting
1. Indentation: 4 spaces
2. Line length: ≤ 80 characters
3. Spacing:
   ```matlab
   % Good
   function result = calculate(a, b)
       result = a + b;
   end
   
   % Bad
   function result=calculate(a,b)
   result=a+b;
   end
   ```

### Comments and Documentation
1. Class documentation:
   ```matlab
   % CLASSNAME Brief description
   %   Detailed description
   %
   % Example:
   %   obj = ClassName(param1, param2);
   %   result = obj.method1(input);
   ```

2. Method documentation:
   ```matlab
   function output = method1(obj, input)
       % METHOD1 Brief description
       %   Detailed description
       %
       % Parameters:
       %   input - Description
       %
       % Returns:
       %   output - Description
       %
       % Example:
       %   output = obj.method1(input);
       
       % Implementation
   end
   ```

## Testing Guidelines

### Test Structure
1. Unit tests:
   ```matlab
   classdef ModalAnalysisTest < matlab.unittest.TestCase
       properties
           TestSurface
       end
       
       methods(TestMethodSetup)
           function setupTest(testCase)
               % Setup code
           end
       end
       
       methods(Test)
           function testFrequencyCalculation(testCase)
               % Test code
           end
       end
   end
   ```

2. Integration tests:
   ```matlab
   function tests = integrationTests
       tests = functiontests(localfunctions);
   end
   
   function testFullAnalysis(testCase)
       % Test code
   end
   ```

### Test Coverage
- Aim for 80% code coverage
- Test edge cases and error conditions
- Include performance tests for critical functions

## Documentation

### Required Documentation
1. Class and function headers
2. Example usage
3. Parameter descriptions
4. Return value descriptions
5. Error conditions
6. Performance considerations

### Documentation Style
```matlab
% FUNCTIONNAME Brief one-line description
%   Detailed description, including mathematical formulation
%   if applicable.
%
% Syntax:
%   output = functionName(input1, input2)
%
% Parameters:
%   input1 - Description [units]
%   input2 - Description [units]
%
% Returns:
%   output - Description [units]
%
% Example:
%   x = 1;
%   y = 2;
%   z = functionName(x, y);
%
% See also: RELATEDFUNCTION1, RELATEDFUNCTION2
%
% References:
%   [1] Author, "Paper Title", Journal, Year
%
% Notes:
%   Additional implementation details or limitations
```

## Pull Request Process

### Before Submitting
1. Run all tests
2. Update documentation
3. Check code style
4. Add examples if needed

### PR Template
```markdown
## Description
Brief description of changes

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Documentation update
- [ ] Performance improvement
- [ ] Code refactoring

## Testing
- [ ] Unit tests added/updated
- [ ] Integration tests added/updated
- [ ] All tests passing

## Documentation
- [ ] Documentation updated
- [ ] Examples added/updated
- [ ] Comments added/updated

## Additional Notes
Any additional information
```

### Review Process
1. Code review by at least one maintainer
2. All tests must pass
3. Documentation must be complete
4. Code style must be consistent

## Questions?
Feel free to open an issue or contact the maintainers.
