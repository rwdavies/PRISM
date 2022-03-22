PRISM (Probabilistic inference of sequence motifs)
==================================================

PRISM is an R package for identifying sequence motifs from hotspots

# Table of contents
1. [Introduction](#paragraph-introduction)
2. [Installation](#paragraph-installation)
3. [Example](#paragraph-example)
4. [Example](#paragraph-tests)

## Introduction <a name="paragraph-introduction"></a>

PRISM is an R package for identifying sequence motifs from hotspots

## Installation <a name="paragraph-installation"></a>

PRISM can be installed in the following way

```
git clone --recursive https://github.com/rwdavies/PRIMSM.git
cd PRISM
./scripts/install-dependencies.sh
./scripts/build-and-install.R
```

## Example <a name="paragraph-example"></a>

An example can be worked through with `example.R` in this directorn

## Tests <a name="paragraph-tests"></a>

Tests are available in `PRISM/tests/testthat/`, and can be run using
```
./scripts/test-unit.sh
./scripts/test-acceptance.sh
```
A single test file of the format `PRISM/tests/testthat/test-<type>-<suffix>.R` can be run using `./scripts/test-<type>.sh <suffix>`, for example `./scripts/test-acceptance.sh one`

See `PRISM/tests/testthat/test-acceptance-one.R` for an example of `getMotifs.R` running to detect motifs using 10,000 hotspots, many with a motif, and 10,000 coldspots without a motif