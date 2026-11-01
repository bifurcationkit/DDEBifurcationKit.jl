using Pkg
cd(@__DIR__)
pkg" activate ."

# pkg"dev DDEBifurcationKit BifurcationKit"
pkg"dev DocumenterCodeBlocks"
# using DocThemeIndigo
ENV["GKSwstype"] = "100"

using Documenter, DDEBifurcationKit, BifurcationKit
using DocumenterCodeBlocks
# to display progress
ENV["JULIA_DEBUG"] = Documenter

makedocs(
	modules = [DDEBifurcationKit, BifurcationKit],
	# pagesonly = true,
	doctest = false,
	sitename = "Bifurcation Analysis of DDEs in Julia",
	warnonly = true,
	draft = false,
	authors = "Romain Veltz",
	plugins = [CodeBlocks()],
	format = Documenter.HTML(
		collapselevel = 1,
		assets=[
			asset("https://bifurcationkit.github.io/assets/js/documentation.js"),
			asset("https://bifurcationkit.github.io/assets/css/documentation.css"),
				],
	),
	pages = Any[
		"🏠 Home" => "index.md",
		"📎 Tutorials" => "tutorials/tutorials.md",
		"🧩 Problems" => [
			"Bifurcation Problem" => "BifProblem.md",
			"Periodic Orbits" => [
				"Introduction" => "periodicOrbit.md",
				"Collocation" => "periodicOrbitCollocation.md",
				],
		],
		"🧰 Functionalities" => [
			"Bifurcations (equilibria)" => [
				"Bifurcation detection (codim 1)" => "detectionBifurcation.md",
				"Fold / Hopf Continuation (codim 2)" => "codim2Continuation.md",
				# "Bogdanov-Takens refinement (codim 3)" => "codim3Continuation.md",
				],
			"Bifurcations (periodic orbits)" => [
				"Bifurcation detection (1 param)" => "detectionBifurcationPO.md",
				# "Fold continuation (2 params)" => "FoldContinuationPO.md",
				# "Period-Doubling continuation (2 params)" => "PDContinuationPO.md",
				# "Neimark-Sacker continuation (2 params)" => "NSContinuationPO.md",
				],
			"Normal form (equilibria)" => [
				"Simple branch point" => "simplebp.md",
				"Non-simple branch point" => "nonsimplebp.md",
				"Simple Hopf" => "simplehopf.md",
			"Normal form (periodic orbit)" => [],
			],
			"Branch switching" => [
				"From equilibria to equilibria" => "abs-from-eq.md",
				"From Hopf/PD/Branch to periodic orbits" => "abs-from-hopf.md",
			],
		],
		"⚙️ Options" => [
			# "Linear Solvers" => "linearsolver.md",
			# "Bordered linear solvers" => "borderedlinearsolver.md",
			"Eigen Solvers" => "eigensolver.md",
			# "Bordered arrays" => "Borderedarrays.md",
		],
		"✏️ Contributing" => [
		],
		"❓ Frequently Asked Questions" => "faq.md",
		"📚 Library" => "library.md"
	]
	)

deploydocs(
	repo = "github.com/bifurcationkit/DDEBifurcationKit.jl.git",
	push_preview = true,
	devbranch = "main"
)
