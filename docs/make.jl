# we use this hacky way because AsymptoticNumericalMethod is not registered
using Pkg

using Documenter, DDEBifurcationKit, Setfield, BifurcationKit
# using DocThemeIndigo
ENV["GKSwstype"] = "100"

# to display progress
# ENV["JULIA_DEBUG"] = Documenter

makedocs(doctest = false,
	sitename = "Bifurcation Analysis of DDEs in Julia",
	format = Documenter.HTML(collapselevel = 1,assets = ["assets/indigo.css"]),
	# format = DocumenterLaTeX.LaTeX(),
	authors = "Romain Veltz",
	pages = Any[
		"Home" => "index.md",
		"Tutorials" => "tutorials/tutorials.md",
		"Problems" => [
			"Bifurcation Problem" => "BifProblem.md",
			"Periodic Orbits" => [
				"Introduction" => "periodicOrbit.md",
				"Collocation" => "periodicOrbitCollocation.md",
				],
		],
		"Functionalities" => [
			"Bifurcations" => [
				"Bifurcation detection (codim 1)" => "detectionBifurcation.md",
				"Fold / Hopf Continuation (codim 2)" => "codim2Continuation.md",
				# "Bogdanov-Takens refinement (codim 3)" => "codim3Continuation.md",
				],
			"Normal form" => [
				"Simple branch point" => "simplebp.md",
				"Non-simple branch point" => "nonsimplebp.md",
				"Simple Hopf" => "simplehopf.md",
			"Normal form (periodic orbit)" => [],
			],
			"Branch switching" => "branchswitching.md",
		],
		"Options" => [
			# "Linear Solvers" => "linearsolver.md",
			# "Bordered linear solvers" => "borderedlinearsolver.md",
			"Eigen Solvers" => "eigensolver.md",
			# "Bordered arrays" => "Borderedarrays.md",
		],
		"Contributing" => [
		],
		"Frequently Asked Questions" => "faq.md",
		"Library" => "library.md"
	]
	)

deploydocs(
	repo = "github.com/bifurcationkit/BifurcationKitDocs.jl.git",
	devbranch = "main"
)
