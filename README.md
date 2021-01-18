# RELAY: a repair algorithm for composite laminates to satisfy lay-up design guidelines


--------------------------------------------------------------------------

The files correspond to the article in the Journal Composites Structures: 
https://doi.org/10.1016/j.compstruct.2020.113448 (2021)

RELAY is a deterministic method for repairing composite laminate lay-ups to satisfy 
stacking-sequence design guidelines.

--------------------------------------------------------------------------
Requirements:

1. A python IDE running python 3.7+

2. The following libraries should accessible:

	- matplotlib
	- pandas
	- numpy

---------------------------------------------------------------------------
In order to use it:

1. clone or download the repository in your desired folder.

2. Set up your environment with the appropriate libraries.

3. Change the settings and run one of the files used to test RELAY: 
run_repair_in_plane.py, run_repair_out_of_plane.py or run_repair.py.

Refer to the documentation for more detail.
--------------------------------------------------------------------------
Folder Structure

- src and subfolders contain the files needed to run the code

- pop contains the files storing the stacking sequence populations used for testing RELAY

- results-paper-RELAY contains the results and analyses generated for the paper 

- run_repair_in_plane.py is used for testing RELAY in-plane steps (steps 1 and 2)

- run_repair_out_of_plane.py is used for testing RELAY out-of-plane steps (steps 3 and 4)

- run_repair.py is used for testing RELAY (steps 1, 2, 3 and 4)

--------------------------------------------------------------------------
Version 1.0.0

--------------------------------------------------------------------------
License

This project is licensed under the MIT License. See the LICENSE for details.

--------------------------------------------------------------------------
Acknowledgments

- Terence Macquart, Paul Weaver and Alberto Pirrera, my PhD supervisors.

--------------------------------------------------------------------------
Author:

Noémie Fedon

For any questions or feedback don't hesitate to contact me at: noemiefedon@gmail.com
or through github at https://github.com/noemiefedon/RELAY
