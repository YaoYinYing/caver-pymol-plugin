CAVER_PLUGIN_VERSION=3.0.3
PROJECT=Caver4


install:
	cp -r Caver4 /Users/yyy/.pymol/startup/

black:
	pre-commit run --all-files

setup-display-gha:
	sudo apt install -y libxkbcommon-x11-0 libxcb-icccm4 libxcb-image0 libxcb-keysyms1 libxcb-randr0 libxcb-render-util0 libxcb-xinerama0 libxcb-xfixes0 x11-utils
	/sbin/start-stop-daemon --start --quiet --pidfile /tmp/custom_xvfb_99.pid --make-pidfile --background --exec /usr/bin/Xvfb -- :99 -screen 0 1920x1200x24 -ac +extension GLX


prepare-test:
	python -m pip install pytest 'pytest-qt>=4.4.0' 'pytest-cov<=6.0.0' 'pytest-mock<3.14.1' pytest-runner 'pytest-xdist>=3.6.1' 'pytest<=8.3.3' pytest-order pytest-emoji pytest-github-actions-annotate-failures 'shellcheck-py==0.10.0.1' pre-commit 'virtualenv<=20.27.1' coverage httpx psutil

test:
	python -m pytest --cov-config=.coveragerc --cov-report=term-missing -v --pyargs --durations=0 -vvv --emoji --cov=$(PROJECT) ./tests


test-faulthandler:
	# see: https://blog.xmatthias.com/post/pytest-debug-segfault/
	python -X faulthandler -m pytest  -p no:faulthandler --cov-config=.coveragerc --cov-report=term-missing -v --pyargs --durations=0 -vvv --emoji --cov=$(PROJECT) ./tests