## R CMD check results

0 errors | 0 warnings | 1 note

The package uses `mirt:::DerivTheta()`, an internal function of the 'mirt'
package. This generates a NOTE about triple-colon access to an unexported
symbol.

`DerivTheta` computes analytic derivatives of item category probabilities with
respect to the latent trait vector. It is used by `population_discrimination()`
to evaluate the gradient of the expected item score for any mirt-supported item
type (2PL, 3PL, GRM, etc.). There is no equivalent public API in 'mirt' that
provides this functionality.

`DerivTheta` has been present and stable in mirt since version 1.0 and is used
by other CRAN packages in the mirt ecosystem (e.g., 'mirtCAT', 'sirt'). The
use is confined to a single call site and is documented in the function's
roxygen comment with a `# nolint: triple_colon_linter` inline suppression to
signal that the use is intentional.
