%% MATLAB Code For Finite Element Analysis
% Problem 4: Determining the stresses and the displacement in a system
% consisting of a 2D truss.
% Units are all in N/mm, mm, or N/mm2!
% Change the code in accordance to the instructions provided in the code
% below.
% Don't forget to change the axes in Function 3:
% plotDeformedTruss.

% Clear all memory from the system
clear;

%% Material Properties - change this code for other trusses.
% E: modulus of elasticity
% A: area of cross section
% L: length of bar
system.E = 70000; 
system.A = 300; 
system.EA = system.E*system.A;

%% Scale factor of the deformed truss - change this code for other trusses
% It's hard to see the deformed shape, so let's scale deformations by a
% factor of 50.
system.scaleFactor = 50;

%% Properties of the structure - change this code for other trusses.
% 0th STEP: Set up a prescribed DOF vector that consists of the DOFs that
% are fixed (zero displacements, i.e. at the supports).
system.prescribedDOFs = [1; 2; 9; 10];

% 1st STEP: Define element node connections, number of nodes and elements.
% Then define the number of DOF's here too.

system.nodeConnections = [1 2; 1 3; 2 3; 2 4; 1 4; 3 4; 3 6; 4 5; 4 6; 3 5; 5 6];
system.numberOfNodes = max(system.nodeConnections,[],'all');
system.numberOfElements = size(system.nodeConnections,1);
system.numberOfDOFs = 2*system.numberOfNodes;

% Define the coordinates of the nodes.
% 1st row: coordinates of node 1: 0,0
% 2nd row: coordinates of node 2: 0,3000
% 3rd row: coordinates of node 3: 3000,0
% 4th row: coordinates of node 4: 3000,3000
% 5th row: coordinates of node 5: 6000,0
% 6th row: coordinates of node 6: 6000,3000
system.coordinatesOfNodes = [0 0; 0 3000; 3000 0; 3000 3000; 6000 0; 6000 3000];

% Define the DOF's for each node. Let's set up a matrix consisting of
% each element's DOF's.
% 1st row: Node 1 DOF's
% 2nd row: Node 2 DOF's
% 3rd row: Node 3 DOF's
% 4th row: Node 4 DOF's
% 5th row: Node 5 DOF's
% 6th row: Node 6 DOF's
system.nodeDOFMatrix = zeros(system.numberOfNodes, 2);
for node = 1:system.numberOfNodes
    system.nodeDOFMatrix(node,:) = [(2*node - 1) (2*node)];
end

%% Solution vectors - change code for the known force vectors for a different truss.
% 2nd STEP: Initiate the load and displacement 
% vectors as zero matrices to be filled in later.
% Load and displacement matrices.
system.displacements = zeros(system.numberOfDOFs, 1);
system.forceVector = zeros(system.numberOfDOFs, 1);

% Each DOF will undergo a displacement and will have a load acting on it.
% The rule that unknown DOF's have known force and DOF = 0 have unknown
% force still applies. Let's include the unknown forces that aren't 0.
system.forceVector(4,1) = -50000;
system.forceVector(8,1) = -100000;
system.forceVector(12,1) = -50000;

%%
% 3rd STEP: Compute the truss stiffness matrix using a function.
system.stiffnessMatrix = formulateTrussStiffnessMatrix(system);

%%
% 4th STEP: Set the boundary conditions.
% The DOF's in system.nodeDOFMatrix that aren't present in
% system.prescribedDOFs are the unknown DOF's.
system.unknownDOFs = setdiff(system.nodeDOFMatrix, system.prescribedDOFs);

%%
% 5th STEP: Determine the displacements, i.e. solve for the unknown DOF's.
% Exact same procedure as in the spring system case.
system.displacements(system.unknownDOFs, 1) = getUnknownDisplacements(system);

%%
% 6th STEP: Plot the truss and its deformed shape.
plotUndeformedTruss(system);
plotDeformedTruss(system);

%%
% 7th STEP: Determine the stresses and forces in each truss element.
system.stressesInElements = getStresses(system);
system.forceVector = getElementForces(system);

%%
% 8th STEP: Display the stresses, forces and displacements.
displayStressesForcesDisplacements(system);

%%
% Function 1: formulateTrussStiffnessMatrix
% This function computes the stiffness matrix of any truss, given its
% structure of arrays is known. The structure's arrays must be renamed as
% follows: nodeConnections, numberOfNodes, numberOfElements, numberOfDOFs,
% coordinatesOfNodes, nodeDOFMatrix, displacements, forceVector (see code
% above line 54).
function [K_Truss_all] = formulateTrussStiffnessMatrix(structure)
% Initiate the truss stiffness matrix as a zeros matrix of size
% structure.numberOfDOFs x structure.numberOfDOFs.
K_Truss_all = zeros(structure.numberOfDOFs, structure.numberOfDOFs);

for element = 1:structure.numberOfElements
    % Determine the cosine, sine and length of each element.
    % Start by determining the coordinates of the elemnt's node
    % coordinates.
    elementNodes = structure.nodeConnections(element,:);
    firstNodeCoordinates = structure.coordinatesOfNodes(elementNodes(1,1), :);
    secondNodeCoordinates = structure.coordinatesOfNodes(elementNodes(1,2), :);
    
    % Determine the length of the element.
    length = sqrt(((firstNodeCoordinates(1,1) - secondNodeCoordinates(1,1))^2)+((firstNodeCoordinates(1,2) - secondNodeCoordinates(1,2))^2));
    
    % Determine cosine and sine.
    c = (secondNodeCoordinates(1,1) - firstNodeCoordinates(1,1))/length;
    s = (secondNodeCoordinates(1,2) - firstNodeCoordinates(1,2))/length;
    
    % Determine the global stiffness matrix of the element.
    K_element_global = ((structure.EA)/length)*[c^2 c*s -c^2 -c*s; c*s s^2 -c*s -s^2;  -c^2 -c*s c^2 c*s; -c*s -s^2 c*s s^2];
    
    % Determine our element's DOF's from our DOF Matrix.
    elementDOFs = [structure.nodeDOFMatrix(elementNodes(1,1),:)'; structure.nodeDOFMatrix(elementNodes(1,2),:)'];
    
    % Determine the all matrix for the element by initiating it as a matrix
    % of zeros of size structure.numberOfNodes x structure.numberOfNodes.
    K_element_all = zeros(structure.numberOfDOFs, structure.numberOfDOFs);
    
    % Then incorporate the global stiffness matrix into the all matrix.
    K_element_all(elementDOFs, elementDOFs) = K_element_global;
    
    % Add the all matrices of all elements together.
    K_Truss_all = K_Truss_all + K_element_all;
end

end

%%
% Function 2: getUnknownDisplacements
% This function takes a truss's structure of arrays and returns the
% displacement of the unknown DOF's.
function [unknownDisplacements] = getUnknownDisplacements(structure)
unknownDisplacements = structure.stiffnessMatrix(structure.unknownDOFs, structure.unknownDOFs)\structure.forceVector(structure.unknownDOFs, 1);
end

%%
% Function 3: plotUndeformedTruss
% This function takes a truss's structure of arrays and plots the
% undeformed truss.
function plotUndeformedTruss(structure)
hold on;
axis([-10, 130, -10, 130]);

for element = 1:structure.numberOfElements
    % Define the nodes of each element
    elementNodes = structure.nodeConnections(element,:);
    firstPoint = structure.coordinatesOfNodes(elementNodes(1,1),:);
    secondPoint = structure.coordinatesOfNodes(elementNodes(1,2),:);
    
    plot([firstPoint(1,1) secondPoint(1,1)],[firstPoint(1,2) secondPoint(1,2)], '-ob', 'LineWidth', 1);
    
end

end

%% Change the axes as required in line 184.
% Function 3: plotDeformedTruss
% This function takes a truss's structure of arrays and plots the deformed
% truss.
function plotDeformedTruss(structure)
hold on;
axis([-10, 6100, -3000, 3100]);

for element = 1:structure.numberOfElements
    % Determine the nodes of each element
    elementNodes = structure.nodeConnections(element,:);
    
    % Determine the DOF's of both nodes
    firstNodedof = structure.nodeDOFMatrix(elementNodes(1,1),:);
    secondNodedof = structure.nodeDOFMatrix(elementNodes(1,2),:);
    
    % Coordinates of each node.
    firstPoint = structure.coordinatesOfNodes(elementNodes(1,1),:)+structure.scaleFactor*[structure.displacements(firstNodedof(1,1),1) structure.displacements(firstNodedof(1,2),1)];
    secondPoint = structure.coordinatesOfNodes(elementNodes(1,2),:)+structure.scaleFactor*[structure.displacements(secondNodedof(1,1),1) structure.displacements(secondNodedof(1,2),1)];
    
    % Plot the deformed shape.
    plot([firstPoint(1,1) secondPoint(1,1)], [firstPoint(1,2) secondPoint(1,2)], '--oc','Linewidth',1);
end
end

%%
% Function 4: getStresses
% This function takes a truss's structure of arrays and returns the
% stresses in each element as a vector.
function [elementStresses] = getStresses(structure)

% Initiate the stresses vector as a zeros vector of size
% structure.numberOfElements x 1.
elementStresses = zeros(structure.numberOfElements,1);

for element = 1:structure.numberOfElements
    % Determine the coordinates of the element's nodes.
    elementNodes = structure.nodeConnections(element,:);
    firstNodeCoordinates = structure.coordinatesOfNodes(elementNodes(1,1),:);
    secondNodeCoordinates = structure.coordinatesOfNodes(elementNodes(1,2),:);
    
    % Determine the length of each element
    length = sqrt(((firstNodeCoordinates(1,1) - secondNodeCoordinates(1,1))^2)+((firstNodeCoordinates(1,2) - secondNodeCoordinates(1,2))^2));
    
    % Determine the cosine and sine of the element.
    c = (secondNodeCoordinates(1,1) - firstNodeCoordinates(1,1))/length;
    s = (secondNodeCoordinates(1,2) - firstNodeCoordinates(1,2))/length;
    
    % Determine the DOF's connected to the element.
    firstNodedof = structure.nodeDOFMatrix(elementNodes(1,1),:);
    secondNodedof = structure.nodeDOFMatrix(elementNodes(1,2),:);
    elementsDOFs = [firstNodedof secondNodedof];
    
    % Determine the stress in each element.
    elementStress = ((structure.E/length)*[-c -s c s])*(structure.displacements(elementsDOFs));
    
    % Input the element's stress into the stress matrix.
    elementStresses(element,1) = elementStress;
end

end

%%
% Function 5: getElementForces
% This function takes a truss's structure of arrays and returns a vector
% consisting of the forces at each node.
function [nodeLoads] = getElementForces(structure)
nodeLoads = structure.stiffnessMatrix*structure.displacements;
end

%%
% Function 6: displayStressesForcesDisplacements
% This function takes a truss's structure of arrays and displays a table
% consisting of the truss's stresses, forces and displacements.
function displayStressesForcesDisplacements(structure)
% Create a row vector of nodes and name it according to the heading on your
% table.
% To name a row vector 'Nodes', we first need to create a row vector
% consisting of nodes from 1 to structure.numberOfDOFs.
% Then we need to convert it into a column vector for it to display in the
% table.
totalNodes = 1:1:structure.numberOfDOFs;
Nodes = totalNodes';

% Rename the force and displacement vectors according to the headings in
% your table.
forces = structure.forceVector;
displacements = structure.displacements;

% Display the table.
T_x = table(Nodes, forces, displacements);
disp('forces and displacements');
disp(T_x);

% Do the exact same thing for the stresses. This time we're working with
% elements instead.
% Again, to display the elements, first we need to create a row vector 
% from 1 to structure.numberOfElements. Then we need to convert it into a
% column vector for it to display in the table.
elementsInOurTruss = 1:1:structure.numberOfElements;
Element = elementsInOurTruss';

Stresses = structure.stressesInElements;
T_y = table(Element,Stresses);
disp('stresses');
disp(T_y);
end









