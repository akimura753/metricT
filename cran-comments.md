## Submission
This is a new submission.

## Test environments
* Ubuntu 24.04, R 4.3.3 (local): R CMD check --as-cran
* Windows, R 4.5: installed from source; the analyses of the accompanying paper were run
* win-builder, R-devel (2026-10-05 r90641), Windows Server 2022: 1 NOTE (see below)

## R CMD check results
0 errors | 0 warnings | 1 note: "New submission", with possibly misspelled words in DESCRIPTION.

## Notes for the reviewer
* The words reported as possibly misspelled (Kimura, Westfall, wPLI, MSC) are
  names of persons or abbreviations defined in the Description text.
* The reference in DESCRIPTION is the published description of the method:
  Kimura (2026) <doi:10.1016/j.neuri.2026.100286>.
* Functions with a `seed` argument call set.seed() only with the value supplied
  through that argument, so that Monte Carlo p-values are reproducible; exact
  (enumerated) tests do not use the random number generator.
