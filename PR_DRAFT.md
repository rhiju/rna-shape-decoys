# Pull Request Drafts - SGNM/DLU/Ciffy Fixes

## Summary
This document contains draft PRs to fix installation issues encountered when setting up SGNM and dependencies for RNA SHAPE prediction from 3D structures.

---

## PR #1: ciffy - Fix duplicate fallback_version

**Repository**: https://github.com/hmblair/ciffy  
**Branch**: main  
**Status**: Ready - commit f7ec301

### Problem
TOML parsing error when installing from source:
```
tomllib.TOMLDecodeError: Cannot overwrite a value (at line 67, column 32)
```

### Root Cause
File: `pyproject.toml` in `[tool.setuptools_scm]` section  
Lines 64 and 67 both had `fallback_version` keys (duplicate keys not allowed in TOML)

### Solution
```diff
[tool.setuptools_scm]
version_scheme = "guess-next-dev"
local_scheme = "node-and-date"
-fallback_version = "0.0.0"
+fallback_version = "0.0.0.dev0"
tag_regex = "^v(?P<version>.*)$"
write_to = "ciffy/_version.py"
-fallback_version = "0.0.0.dev0"
```

### PR Details
- **Title**: Fix: Remove duplicate fallback_version in setuptools_scm configuration
- **Type**: Bug fix
- **Affects**: All users installing ciffy from source
- **Tested**: ✅ ciffy imports successfully after fix

---

## PR #2: sgnm - Update dlu dependency and fix imports

**Repository**: https://github.com/hmblair/sgnm  
**Branch**: master  
**Status**: Ready - commit pending

### Problem #1: Incorrect dlu import
File: `sgnm/models.py` line 14

**Current** (broken):
```python
from dlu import RadialBasisFunctions, DenseNetwork
```

**Should be**:
```python
from dlu.modules import RadialBasisFunctions, DenseNetwork
```

**Reason**: dlu uses namespace package structure; classes are in `dlu.modules` submodule

### Problem #2: Dependency name mismatch
File: `pyproject.toml` line 20

**Current**:
```toml
dependencies = [
    ...
    "dlu",
]
```

**Should be**:
```toml
dependencies = [
    ...
    "dlu-torch",
]
```

**Reason**: Package published as `dlu-torch` on PyPI (hmblair/dlu repository)

### PR Details
- **Title**: Fix: Update dlu import and dependency name
- **Type**: Bug fix / Maintenance  
- **Affects**: All users installing sgnm from source
- **Changes**:
  1. Update models.py import: `from dlu import` → `from dlu.modules import`
  2. Update pyproject.toml: `"dlu"` → `"dlu-torch"`
- **Tested**: ✅ SGNM models import successfully after fixes

### Detailed Changes

```diff
# sgnm/models.py
-from dlu import RadialBasisFunctions, DenseNetwork
+from dlu.modules import RadialBasisFunctions, DenseNetwork

# pyproject.toml
dependencies = [
    "torch>=2.0",
    "ciffy",
    "h5py",
    "tqdm",
-   "dlu",
+   "dlu-torch",
]
```

---

## PR #3: dlu - Documentation update (optional enhancement)

**Repository**: https://github.com/hmblair/dlu  
**Status**: Recommended - not critical

### Suggestion
Add note to README clarifying:
1. Package is published as `dlu-torch` on PyPI
2. Imported as `dlu` in code
3. Classes are in `dlu.modules` submodule (for namespace compatibility)

Example README note:
```markdown
## Installation

Install via pip:
```bash
pip install dlu-torch
```

Import in code:
```python
from dlu import DenseNetwork, RadialBasisFunctions
# or explicitly:
from dlu.modules import DenseNetwork, RadialBasisFunctions
```
```

---

## Additional Issue: Model Checkpoint Distribution

**Problem**: Automated downloads of SGNM model checkpoints fail (return 301 redirects)

**Current State**: 
- Users cannot easily download model weights via `wget` or `curl`
- Checkpoints hosted on GitHub releases but redirect URLs don't work

**Recommendation for hmblair/sgnm**:
1. Add GitHub CLI download instructions to README:
   ```bash
   gh release download v1.0.0 -p "gnm-checkpoint.pth" -D .
   ```

2. Or provide direct HuggingFace Hub URLs if permissions allow

3. Document manual download location clearly

---

## Implementation Notes

All changes have been tested in the context of using SGNM for RNA SHAPE prediction:
- ✅ ciffy: Installs and imports successfully
- ✅ dlu-torch: Exports all required classes
- ✅ sgnm: Models.py imports and loads successfully
- ✅ Prediction scripts: Ready to run with model checkpoints

**Repo Status**: All repos cloned from correct hmblair sources and on master/main branches
- sgnm: https://github.com/hmblair/sgnm.git (master)
- dlu: https://github.com/hmblair/dlu.git (main)  
- ciffy: https://github.com/hmblair/ciffy.git (main)

