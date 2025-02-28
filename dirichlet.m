%{
	Author: Rodrigo Kiefe

	Solves the magnetic field H using the magnetic field scalar 
	potential and the finite element method

	A dirichelet boundary condition (u=0 at the outside) is used to impose
	a unique solution

	Requires the pde tool box from matlab to generate the model and mesh.
	Everything else can be done in whatever programming language you want.

	--- Variables ---
	Hmax -> max element size
	Hmin -> min elemnet size
	Hrefined -> local refinement on the object's surface
%}

close all
clear
clc

% ------------ User defined variables and constants ------------
hmax = 2; 		% 0 -> let matlab choose
hmin = 0.05; 		% ...
Hrefined = 0.05;	% Local refinement
mu0 = pi*4e-7;
% --------------------------------------------------------------

% PDE model
model = geometry(1); % 1 for circle, 2 for rectangle

% >> Mesh
generateMesh(model, ...
             "GeometricOrder","linear", ...
             Hmax=hmax, ...
             Hmin=hmin, ...
             Hgrad=1.5, ...
             Hedge={5:model.Geometry.NumEdges,Hrefined});

mesh = processMesh(model,Hmax=hmax,Hmin=hmin); % 0 for 'let matlab choose'

fprintf("\n n nodes %d\n",mesh.nv)
fprintf("\n n triangles %d\n",mesh.nt)
fprintf("\n nInside %d\n",mesh.nInside)
pause()


% Magnetization
M = zeros(mesh.nInside,2) + [1,0];

% >> Stiffness matrix
A = zeros(mesh.nv,mesh.nv);
for k = 1:mesh.nt
	nds = mesh.t(1:3,k);

	b = zeros(3,1); c = zeros(3,1);
	for i = 1:3
		[~,b(i),c(i)] = abc(mesh.p,nds,nds(i));
	end

	A(nds,nds) = A(nds,nds) + mesh.VE(k)*(b*b' + c*c');
end


% Load Vector
q = zeros(mesh.nv,1);
for ik = 1:mesh.nInside
	k = mesh.InsideElements(ik);
	nds = mesh.t(1:3,k);

	for i = 1:3
		[~,b,c] = abc(mesh.p,nds,nds(i));
		q(nds(i)) = q(nds(i)) + mesh.VE(k)*dot(M(ik,:),[b,c]);
	end
end

% >> Magnetic scalar potential

% Dirichlet boundary condition
fixed = findNodes(model.Mesh,"region","Edge",1:4);
uD = zeros(numel(fixed),1); % Potential is zero at external boundary

% Solve in the free nodes
free = setdiff(1:mesh.nv,fixed);

u = zeros(mesh.nv,1);
u(fixed) = uD;
u(free) = A(free,free)\(q(free)-A(free,fixed)*uD);

pdegplot(model); hold on
scatter(mesh.p(1,:),mesh.p(2,:),[],u,'filled')

% Magnetic field H
H = gradU(mesh,u);

% >> Norm of H
H_mod = sqrt(sum(H.^2,2)); % |H|

% >> Element centroids
pc = zeros(mesh.nt,2);
for k = 1:mesh.nt
	pc(k,:) = mean(mesh.p(1:2,mesh.t(1:3,k)),2);
end

% >> Plots
figure
pdegplot(model); hold on
scatter(pc(mesh.InsideElements,1),...
		pc(mesh.InsideElements,2),[],...
		H_mod(mesh.InsideElements),'filled')

cbar = colorbar;
cbar.Label.String = "|H|";

figure
pdegplot(model); hold on
quiver(pc(mesh.InsideElements,1),...
	   pc(mesh.InsideElements,2),...
	   H(mesh.InsideElements,1),...
	   H(mesh.InsideElements,2))



% Functions
function mesh = processMesh(model,options)

	arguments
		model;
		options.Hmax = 0;
		options.Hmin = 0;
	end
	mesh = struct;

	[~,mesh.VE] = area(model.Mesh); % area of each element

	% Process mesh to get a list of edges
	[mesh.p,mesh.edgeList,mesh.t] = meshToPet(model.Mesh);

	% edge list: nd1, nd2, edge index
	mesh.edgeList = mesh.edgeList([1,2,5],:);

	mesh.nv = length(mesh.p); % Number of nodes
	mesh.nt = length(mesh.t); % number of elements
	mesh.ne = length(mesh.edgeList); % number of edges

	mesh.p = [mesh.p;zeros(1,mesh.nv)];
	for e = 1:mesh.ne
		mesh.p(end,mesh.edgeList(1:2,e)) = mesh.edgeList(end,e);
	end

	mesh.InsideElements = findElements(model.Mesh,"region",Face=2);
	mesh.nInside = numel(mesh.InsideElements);
end

function H = gradU(mesh,u)
	H = zeros(mesh.nt,2);
	for k = 1:mesh.nt
		% Nodes of the element
		nds = mesh.t(1:3,k);
	 
	    % Calculate the total potential of the current element
	    for ind = 1:length(nds)
	        nd = nds(ind);

	        [~,bi,ci] = abc(mesh.p,nds,nd);

	        H(k,1) = H(k,1) - bi*u(nd);
	        H(k,2) = H(k,2) - ci*u(nd);
	    end
	end
end

function b = loadVector(mesh,f)
	b = zeros(mesh.nv,1);
	for k = 1:mesh.nt
		nds = mesh.t(1:3,k);
		for i = 1:length(nds)
			nd = nds(i);
			[~,bi,ci] = abc(mesh.p,nds,nd);

			b(nd) = b(nd) + mesh.AE(k)*dot([bi,ci],f(k,:));
		end
	end
end

function C = Cmaker(mesh)
	C = zeros(1,mesh.nv);
	for k = 1:mesh.nt
	    nds = mesh.t(1:3,k);

	    % For each node of the element
	    for ind = 1:length(nds)
	        nd = nds(ind);

	        % a b c d parameters for that element and node
	        [a,b,c] = abc(mesh.p,nds,nd);

	        % Corrdinate of element center
	        pcenter = mean(mesh.p(1:2,nds),2);

	        C(nd) = C(nd) + (a + b*pcenter(1) + c*pcenter(2))*mesh.AE(k);
	    end
	end
end

function A = stiffnessMatrix(mesh,mu)
	A = zeros(mesh.nv,mesh.nv);
	for k = 1:mesh.nt
		
		% Nodes of the element
		nds = mesh.t(1:3,k);

		% For each node
		for i = 1:length(nds)
			[~,bi,ci] = abc(mesh.p,nds,nds(i)); % basis function parameters

			for j = i:length(nds) % matrix is symetric
				[~,bj,cj] = abc(mesh.p,nds,nds(j));

				A(nds(j),nds(i)) = A(nds(j),nds(i)) + mu(k)*(bi*bj + ci*cj)*mesh.AE(k);
				A(nds(i),nds(j)) = A(nds(j),nds(i)); % A is symetric
			end
		end
	end
end

function [a,b,c] = abc(p,nds,nd) % [a,b,c]
	nd1 = 0;
	nd2 = 0;
	% Find opposing nodes to nd
	for i = 1:length(nds)
		if nds(i) ~= nd && nd1 == 0
			nd1 = nds(i);
		end

		if nd1 ~= nds(i) && nds(i) ~= nd && nd2 == 0
			nd2 = nds(i);
		end
	end
	
	% Vectors in-plane
	v1 = [p(1,nd)-p(1,nd1), p(2,nd)-p(2,nd1),1];
	v2 = [p(1,nd)-p(1,nd2), p(2,nd)-p(2,nd2),1];

	n = cross(v1,v2);

	a = 1 + (n(1)*p(1,nd) + n(2)*p(2,nd))/n(3);
	b = -n(1)/n(3);
	c = -n(2)/n(3);
end