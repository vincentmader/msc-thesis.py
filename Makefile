setup:
	cd bin && ./create_virtualenv.sh
	mkdir -p lib
	cd bin && ./clone_matplotlib_themes.sh
	cd bin && ./clone_radau_solver.sh
	mkdir -p out
