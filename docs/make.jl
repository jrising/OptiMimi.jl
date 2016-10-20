using Documenter, OptiMimi

makedocs(
         modules=[OptiMimi],
         pages = Any[
                     "Introduction" => "index.md",
                     "Linear Programming" => "linproghouse.md",
                     ]
         )
