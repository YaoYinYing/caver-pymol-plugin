from packaging import version

from Caver4.caver_pymol import VERSION


class TestVersionNumber:
    def test_current_version(self):
        version.parse(VERSION)
