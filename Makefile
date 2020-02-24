shared: setup.py barneslib.cc
	mkdir shared
	python setup.py build
	cp */*/*.so shared/barneslib.so
clean:
	rm -rf shared
	rm -rf build
