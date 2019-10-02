# LakeMaggiore

Simulation model of Lake Maggiore that can be combined with the [Borg Multi-Objective Evolutionary Algorithm](http://borgmoea.org/) or another MOEA implemented in the [MOEA Framework](http://moeaframework.org/) to design a set of Pareto optimal operating policies via [Evolutionary Multi-Objective Direct Policy Search](https://ascelibrary.org/doi/abs/10.1061/(ASCE)WR.1943-5452.0000570). 

Lake Maggiore is an subalpine lake whose 6,598 km2 wide watershed is equally divided between Italy and Switzerland. The hydrologic regime in the catchment is typical of Alpine zones with snow-melt driving most of the inflow in spring and early summer and autumn rains producing a second peak around October. The Lake has been regulated since 1943 with the primary purpose of downstream water allocation. The primary water uses include a large irrigated area downstream as well as hydropower production. A seasonal inflow pattern is observed, with the lake often filling twice a year, in late spring and in autumn. Each period is followed by a draw down cycle. This release strategy can conflict with the secondary goal of using the lake storage for buffering flood peaks. Several flood episodes affecting the cities around the lake’s shoreline (about 168.000 inhabitants) have been registered in the last decades, with moderate to severe economic losses in terms of property and activities. At the same time, downstream users have experienced frequent summer water deficits which have consequences for both farmers and hydropower producers, as well as the environment.

The two primary objectives (both to be minimized) are hence formulated as follows:
* Flood control: the average annual number of flooding days in the simulation horizon.
* Water supply deficit: the daily average quadratic water deficit between lake releases and downstream water demands (comprising both farmers and hydropower producers)

To compile and run:
* Run `make` in the test folder to compile
* Run `./LakeMaggioreSim settings_lakeMaggiore.txt < u_test29.txt` to perform a simulation with a random policy


----
### Copyright:
  
Copyright 2019 [NRM group - Politecnico di Milano](www.nrm.deib.polimi.it)
  
Developers: Simona Denaro, Matteo Giuliani, Andrea Castelletti.
  
LakeMaggiore is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
  
The code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
  
You should have received a copy of the GNU General Public License along with LakeMaggiore.  If not, see <http://www.gnu.org/licenses/licenses.en.html>.
