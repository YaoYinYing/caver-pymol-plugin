CAVER_PLUGIN_VERSION=3.0.3

dist:
	mkdir dist
	zip -r dist/caver_${CAVER_PLUGIN_VERSION}.zip README.md LICENSE CHANGELOG COPYING Caver3

clean:
	rm -rf dist

install:
	cp -r Caver4 /Users/yyy/.pymol/startup/

black:
	pre-commit run --all-files