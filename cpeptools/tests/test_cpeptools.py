"""
Unit and regression test for the cpeptools package.
"""

# Import package, test suite, and other packages as needed
import cpeptools
import pytest
import sys

def test_cpeptools_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "cpeptools" in sys.modules
