import setuptools

setuptools.setup(
	name = "uwstatus",
	version = "0.1",
	author = "Ava Polzin",
	author_email = "apolzin@uchicago.edu",
	description = "Match object/coordinates to Dragonfly UW field + see observation status.",
	packages = ["uwstatus", "uwstatus/status", "uwstatus/plotting"],
	url = "https://github.com/avapolzin/UWSStatusSearch",
	license = 'MIT',
	classifiers = [
		"Development Status :: 4 - Beta",
		"Intended Audience :: Science/Research",
		"License :: OSI Approved :: MIT License",
		"Operating System :: OS Independent",
		"Programming Language :: Python",
		"Topic :: Scientific/Engineering :: Astronomy",
		"Topic :: Scientific/Engineering :: Physics"],
	python_requires = ">=3",
	install_requires = ["numpy", "matplotlib", "astropy", "pandas", 
		"datetime", "pygithub", "descartes"]
)
