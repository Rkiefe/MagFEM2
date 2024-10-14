%{
	Magnetostatics with FEM
	
	: - div H = div M
	-> nabla^2 u = div M
	-> int_omega (1+chi) grad u dot grad phi dV = int_omga chi Hext dV  

	Also uses the lagrange multiplier technique to get a unique solution.
%}

clear
close all
clc

mu0 = pi*4e-7;
Happ = 1/mu0 *[0,1];

model = geometry(2);

% % PDE model
% model = createpde();

% % Add a rectangle
% container_width = 2;
% container_height = 1;
% position = [0,0];

% box_size = [container_width,container_height]; % width and hight
% rect = [3;
% 	4;
% 	position(1)-box_size(1)/2;
% 	position(1)+box_size(1)/2;
% 	position(1)+box_size(1)/2;
% 	position(1)-box_size(1)/2;
% 	position(2)-box_size(2)/2;
% 	position(2)-box_size(2)/2;
% 	position(2)+box_size(2)/2;
% 	position(2)+box_size(2)/2];
% gd = rect;

% % Names for the shapes
% ns = char('rect1');
% ns = ns';

% % Set logic for combining shapes
% sf = 'rect1'; 

% % Create geometry
% [dl,bt] = decsg(gd,sf,ns);
% geometryFromEdges(model,dl);

% >> Plot model
pdegplot(model); hold on

% >> Mesh
mesh = processMesh(model);

% Hext
Hext = zeros(mesh.nt,2) + Happ;		% applied field
chi = zeros(mesh.nt,1);					% susceptibility
chi(mesh.InsideElements) = 3;

% >> Stiffness matrix
A = stiffnessMatrix(mesh,1+chi);

% Lagrange multiplier tec.
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

% Load Vector
b = loadVector(mesh,chi.*Hext);


mat = zeros(mesh.nv+1);
mat(1:mesh.nv,1:mesh.nv) = A;
mat(end,1:mesh.nv) = C;
mat(1:mesh.nv,end) = C';

u = mat\[b;0];
u = u(1:end-1);

% H
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

H = H + Happ;

% >> Norm of H
H_mod = sqrt(sum(H.^2,2)); % |H|



pc = zeros(mesh.nt,2);
for k = 1:mesh.nt
	pc(k,:) = mean(mesh.p(1:2,mesh.t(1:3,k)),2);
end

scatter(pc(:,1),pc(:,2),10,mu0*H_mod,'filled')



% Functions
function mesh = processMesh(model,options)

	arguments
		model;
		options.Hmax = 0;
		options.Hmin = 0;
	end

	generateMesh(model,"GeometricOrder","linear","Hmax",options.Hmax,"Hmin",options.Hmin);

	[~,AE] = area(model.Mesh); % area of each element

	% Process mesh to get a list of edges
	[p,edgeList,t] = meshToPet(model.Mesh);

	% edge list: nd1, nd2, edge index
	edgeList = edgeList([1,2,5],:);

	nv = length(p); % Number of nodes
	nt = length(t); % number of elements
	nedg = length(edgeList); % number of edges

	p = [p;zeros(1,nv)];
	for e = 1:nedg
		p(end,edgeList(1:2,e)) = edgeList(end,e);
	end

	
	mesh = struct;
	mesh.p = p; clear p
	mesh.t = t; clear t
	mesh.nv = nv; clear nv
	mesh.nt = nt; clear nt
	mesh.ne = nedg; clear nedg
	mesh.edgeList = edgeList; clear edgeList
	mesh.AE = AE; clear AE
	mesh.InsideElements = findElements(model.Mesh,"region",Face=2);
	mesh.nInside = numel(mesh.InsideElements);

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