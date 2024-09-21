local Pathfinding = {}

local Workspace = game:GetService("Workspace")
local RunService = game:GetService("RunService")

local mainScript = script.Parent
local utilFunctions = require(mainScript:WaitForChild("UtilFunctions"))
local voxelMap = require(mainScript:WaitForChild("VoxelMap"))
local priorityQueue = require(script:WaitForChild("PriorityQueue"))

local FACE_VERTICES = 4
local SHIP_SPEED = 15
local NODE_SIZE = Vector3.new(0.2, 0.2, 0.2)
local MIN = math.min
local DIMENSIONS = {"X", "Y", "Z"}

local allLengths = {}

local function nodeToString(node)
	return "("..node.Index[1]..", "..node.Index[2]..", "..node.Index[3]..")"
end

--returns true if first set of keys is greater then second set of keys, key = {length + heuristics, length}
local comparator = function(priorityOne, priorityTwo)
	if priorityOne[1] > priorityTwo[1] or (priorityOne[1] == priorityTwo[1] and priorityOne[2] > priorityTwo[2]) then
		return true
	else
		return false
	end
end

--get length type of a node for the given ship (LookAhead, Goal)
local function getLength(ship, node, lengthType)
	local shipLengths = allLengths[ship]
	local nodeLengths = shipLengths and shipLengths[node]
	if nodeLengths then
		return nodeLengths[lengthType]
	elseif shipLengths and lengthType == "Goal" then
		shipLengths[node] = {Goal = math.huge}
		return math.huge
	end
end

--set length type of a node for the given ship (LookAhead, Goal)
local function setLength(ship, node, lengthType, value)
	local shipLengths = allLengths[ship]
	if shipLengths then
		local nodeLengths = shipLengths[node]
		if not nodeLengths then
			shipLengths[node] = {}
			nodeLengths = shipLengths[node]
		end
		nodeLengths[lengthType] = value
	end
end

--Euclidean distance
local function heuristic(startPosition, endPosition)
	return (startPosition - endPosition).magnitude
end

local function computeKeys(ship, startNode, node)
	local goalLength = getLength(ship, node, "Goal")
	local lookAheadLength = getLength(ship, node, "LookAhead")
	local minimumLength = MIN(goalLength, lookAheadLength)
	local keys = {minimumLength + heuristic(startNode.Position, node.Position), minimumLength}
	
	return keys
end

local function getNearestNode(position)
	
end

local function get1DNeighbors(nodeMap, node)
	local neighbors = {}
	for i = -1, 1, 2 do
		local yList = nodeMap[node.Index[1]]
		local zList = yList and yList[node.Index[2]]
		local adjacentNode = zList and zList[node.Index[3] + i]
		if adjacentNode then
			table.insert(neighbors, adjacentNode)
		end
	end
	return neighbors
end

local function get2DNeighbors(nodeMap, node)
	local neighbors = {}
	if not node then
		warn("Cannot get 2D neighbors - node does not exist")
		return neighbors
	end
	for i = -1, 1, 2 do
		local yList = nodeMap[node.Index[1]]
		local zList = yList and yList[node.Index[2] + i]
		local adjacentNode = zList and zList[node.Index[3]]
		if adjacentNode then
			table.insert(neighbors, adjacentNode)
			local oneDimensionNeighbors = get1DNeighbors(nodeMap, adjacentNode)
			for index, neighbor in pairs(oneDimensionNeighbors) do
				table.insert(neighbors, neighbor)
			end
		end
	end
	local oneDimensionNeighbors = get1DNeighbors(nodeMap, node)
	for index, neighbor in pairs(oneDimensionNeighbors) do
		table.insert(neighbors, neighbor)
	end
	return neighbors
end

local function getNeighboringNodes(nodeMap, node)
	local neighbors = {}
	if not node then
		warn("Cannot get neighbors - node does not exist")
		return neighbors
	end
	for i = -1, 1, 2 do
		local yList = nodeMap[node.Index[1] + i]
		local zList = yList and yList[node.Index[2]]
		local adjacentNode = zList and zList[node.Index[3]]
		if adjacentNode then
			table.insert(neighbors, adjacentNode)
			local twoDimensionNeighbors = get2DNeighbors(nodeMap, adjacentNode)
			for index, neighbor in pairs(twoDimensionNeighbors) do
				table.insert(neighbors, neighbor)
			end
		end
	end
	local twoDimensionNeighbors = get2DNeighbors(nodeMap, node)
	for index, neighbor in pairs(twoDimensionNeighbors) do
		table.insert(neighbors, neighbor)
	end
	return neighbors
end

--multiplying by the index changes the sign appropriately
local function getXFaces(nodeMap, node, index)
	local faces = {}
	local xOrder = {{0, 1, 0}, {0, 0, 1 * index}, {0, -1, 0}, {0, 0, -1 * index}}
	for faceIndex = 1, 4 do
		local face = {}
		table.insert(face, node)
		local currentNode = node
		for directionIndex = 1, 3 do
			local direction = xOrder[directionIndex]
			local yList = nodeMap[currentNode.Index[1] + direction[1]]
			local zList = yList and yList[currentNode.Index[2] + direction[2]]
			local adjacentNode = zList and zList[currentNode.Index[3] + direction[3]]
			if adjacentNode then
				table.insert(face, adjacentNode)
				currentNode = adjacentNode
			else
				break
			end
		end
		local firstDirection = xOrder[1]
		table.remove(xOrder, 1)
		table.insert(xOrder, firstDirection)
		if #face == 4 then
			table.insert(faces, face)
		end
	end
	return faces
end

local function getYFaces(nodeMap, node, index)
	local faces = {}
	local yOrder = {{1, 0, 0}, {0, 0, -1 * index}, {-1, 0, 0}, {0, 0, 1 * index}}
	for i = 1, 4 do
		local face = {}
		table.insert(face, node)
		local currentNode = node
		for i = 1, 3 do
			local direction = yOrder[i]
			local yList = nodeMap[currentNode.Index[1] + direction[1]]
			local zList = yList and yList[currentNode.Index[2] + direction[2]]
			local adjacentNode = zList and zList[currentNode.Index[3] + direction[3]]
			if adjacentNode then
				table.insert(face, adjacentNode)
				currentNode = adjacentNode
			else
				break
			end
		end
		local firstDirection = yOrder[1]
		table.remove(yOrder, 1)
		table.insert(yOrder, firstDirection)
		if #face == 4 then
			table.insert(faces, face)
		end
	end
	return faces
end

local function getZFaces(nodeMap, node, index)
	local faces = {}
	local zOrder = {{0, 1, 0}, {-1 * index, 0, 0}, {0, -1, 0}, {1 * index, 0, 0}}
	for i = 1, 4 do
		local face = {}
		table.insert(face, node)
		local currentNode = node
		for i = 1, 3 do
			local direction = zOrder[i]
			local yList = nodeMap[currentNode.Index[1] + direction[1]]
			local zList = yList and yList[currentNode.Index[2] + direction[2]]
			local adjacentNode = zList and zList[currentNode.Index[3] + direction[3]]
			if adjacentNode then
				table.insert(face, adjacentNode)
				currentNode = adjacentNode
			else
				break
			end
		end
		local firstDirection = zOrder[1]
		table.remove(zOrder, 1)
		table.insert(zOrder, firstDirection)
		if #face == 4 then
			table.insert(faces, face)
		end
	end
	return faces
end

local function getAdjacentFaces(nodeMap, node)
	local faces = {}
	for i = -1, 1, 2 do
		local xYList = nodeMap[node.Index[1] + i]
		local xZList = xYList and xYList[node.Index[2]]
		local adjacentXNode = xZList and xZList[node.Index[3]]
		if adjacentXNode then
			local xFaces = getXFaces(nodeMap, adjacentXNode, i)
			for index, face in pairs(xFaces) do
				table.insert(faces, face)
			end
		end
		
		local yYList = nodeMap[node.Index[1]]
		local yZList = yYList and yYList[node.Index[2] + i]
		local adjacentYNode = yZList and yZList[node.Index[3]]
		if adjacentYNode then
			local yFaces = getYFaces(nodeMap, adjacentYNode, i)
			for index, face in pairs(yFaces) do
				table.insert(faces, face)
			end
		end
		
		local zYList = nodeMap[node.Index[1]]
		local zZList = zYList and zYList[node.Index[2]]
		local adjacentZNode = zZList and zZList[node.Index[3] + i]
		if adjacentZNode then
			local zFaces = getZFaces(nodeMap, adjacentZNode, i)
			for index, face in pairs(zFaces) do
				table.insert(faces, face)
			end
		end
	end
	return faces
end

--[[local function getTraversalCost(voxelMap, startNode, edge)
	local voxelOffset = {0, 0, 0}
	for i = 1, 3 do
		--local offset = startNode.Index[1] - endNode.Index[1]
		--if offset < 0 then
			voxelOffset[i] = -1
		--end
	end
	local yList = voxelMap[startNode.Index[1] + voxelOffset[1]]
	--local zList = yList and yList[startNode.Index[2] + voxelOffset[2]]
	--local adjacentVoxel = zList and zList[startNode.Index[3] + voxelOffset[3]]
	--if adjacentVoxel then
	--	return adjacentVoxel.TraversalCost
	--end
	--return 1
--end
--c is traversal cost of cell with corners s, s1, s2;
--b is traversal cost of cell with corners s, s1 but not s2;
local function sharedDimensions(nodeOne, nodeTwo)
	local sharedCount = 0
	
	local positionOne = nodeOne.Position
	local positionTwo = nodeTwo.Position
	for index, dimension in pairs(DIMENSIONS) do
		if positionOne[dimension] == positionTwo[dimension] then
			sharedCount = sharedCount + 1
		end
	end
	
	return sharedCount
end

local function computeEdgeCost(ship, node, edge)
	local sideOne
	local sideTwo
	--[[local difference = node.Position - edge[1].Position
	local isDiagonal = math.abs(difference.X) > 0 and math.abs(difference.Y) > 0 and math.abs(difference.Z) > 0
	if isDiagonal then
		sideOne = edge[2]
		sideTwo = edge[1]
	else
		sideOne = edge[1]
		sideTwo = edge[2]
	end]]
	if sharedDimensions(node, edge[1]) > sharedDimensions(node, edge[2]) then
		sideOne = edge[1]
		sideTwo = edge[2]
	else
		sideOne = edge[2]
		sideTwo = edge[1]
	end
	local traversalOne = 1
	local traversalTwo = 1
	local sideOneGoal = getLength(ship, sideOne, "Goal")
	local sideTwoGoal = getLength(ship, sideTwo, "Goal")
	local totalCost
	if MIN(traversalOne, traversalTwo) == math.huge then
		totalCost = math.huge
	elseif sideOneGoal <= sideTwoGoal then
		totalCost = MIN(traversalOne, traversalTwo) + getLength(ship, sideOne, "Goal")
	else
		local costDifference = sideOneGoal - sideTwoGoal
		if costDifference <= traversalTwo then
			if traversalOne <= costDifference then
				totalCost = traversalOne * math.sqrt(2) + sideTwoGoal
			else
				local y = MIN(costDifference/math.sqrt((traversalOne * traversalOne) - (costDifference * costDifference)), 1)
				totalCost = (traversalOne * math.sqrt(1 + y * y)) + (costDifference * (1 - y)) + sideTwoGoal
			end
		else
			if traversalOne <= traversalTwo then
				totalCost = traversalOne * math.sqrt(2) + sideTwoGoal
			else
				local x = 1 - MIN(traversalTwo/math.sqrt((traversalOne * traversalOne) - (traversalTwo * traversalTwo)), 1)
				totalCost = (traversalOne * math.sqrt(1 + (1 - x) * (1 - x))) + (traversalTwo * x) + sideTwoGoal
			end
		end
	end
	return totalCost
end

local function computeEdgeMinima(ship, node, face)--face = {s1, s2, s3, s0}
	local minima = {}
	local debuging = false
	
	if node.Index[1] == 1 and node.Index[2] == 3 and node.Index[3] == 2 then
		debuging = true
	end
	if debuging then
		print("Compute Edge Minima for face: "..nodeToString(face[1])..", "..nodeToString(face[2])..", "..nodeToString(face[3])..", "..nodeToString(face[4]))
	end
	for i = 1, FACE_VERTICES do
		local correspondingIndex = i == FACE_VERTICES and 1 or i + 1
		local nodeOne = face[i]
		local nodeTwo = face[correspondingIndex]
		--[[
		local nodeDistance = heuristic(nodeOne.Position, nodeTwo.Position)
		local nodeOneLength = getLength(ship, nodeOne, "Goal")
		local nodeTwoLength = getLength(ship, nodeTwo, "Goal")
		local edgeMinima = nodeDistance * nodeTwoLength + (1 - nodeDistance) * nodeOneLength--pretty sure the '1' is supposed to be voxel size
		--]]
		local edgeMinima = computeEdgeCost(ship, node, {nodeOne, nodeTwo})
		if debuging then
			print("MinimaCost: "..tostring(edgeMinima).." for edge: "..nodeToString(nodeOne)..", "..nodeToString(nodeTwo))
		end
		table.insert(minima, edgeMinima)
	end
	local minimaPrint = ""
	if debuging then
		for i, edgeMinima in pairs(minima) do
			minimaPrint = minimaPrint .. edgeMinima .. " "
		end
	end
	print(minimaPrint)
	return minima
end
--local u0, t1, u1, t0 = unpack(computeEdgeMinima(ship, node, face))--{u0, t1, u1, t0}
local function getTraversalCost(voxelMap, startNode, endNode)
	local voxelOffset = {0, 0, 0}
	for i = 1, 3 do
		local offset = startNode.Index[1] - endNode.Index[1]
		if offset < 0 then
			voxelOffset[i] = -1
		end
	end
	local yList = voxelMap[startNode.Index[1] + voxelOffset[1]]
	local zList = yList and yList[startNode.Index[2] + voxelOffset[2]]
	local adjacentVoxel = zList and zList[startNode.Index[3] + voxelOffset[3]]
	if adjacentVoxel then
		return adjacentVoxel.TraversalCost
	else
		--warn("could not get traversal cost between nodes: "..nodeToString(startNode)..", "..nodeToString(endNode))
	end
	return 1
end

local function computeCost(ship, voxelMap, node, face)
	local u0, t1, u1, t0 = unpack(computeEdgeMinima(ship, node, face))--{u0, t1, u1, t0}
	local xMinima = ((t1 - t0) * u0 + t0) / (1 - (t1 - t0) * (u1 - u0))
	local yMinima = (u1 - u0) * xMinima + u0
	local traversalCost = getTraversalCost(voxelMap, node, face[3])--face[1] will always be the node on the face that shares 2 axis with the origin node, face[3] will be the node on the face that does not share any axis
	local s1Length = getLength(ship, face[1], "Goal")
	local s2Length = getLength(ship, face[2], "Goal")
	local s3Length = getLength(ship, face[3], "Goal")
	local s0Length = getLength(ship, face[4], "Goal")
	local totalCost = traversalCost * math.sqrt(1 + (xMinima * xMinima) + (yMinima * yMinima))
					+ ((s1Length + (s0Length - s1Length)) * xMinima) * (1 - yMinima)
					+ ((s2Length + (s3Length - s2Length)) * xMinima) * yMinima
	return totalCost
end

local function updateState(ship, nodeMap, voxelMap, queue, goalNode, startNode, node, goalAdjacent)
	print("---Update Node: ("..node.Index[1]..", "..node.Index[2]..", "..node.Index[3]..")---")
	local goalLength = getLength(ship, node, "Goal")
	if not goalLength then--have not visited node yet
		setLength(ship, node, "Goal", math.huge)
		goalLength = math.huge
	end
	if node ~= goalNode then--not utilFunctions:InList(goalFace, node) then
		--node.ChangeColor(BrickColor.new("Bright green"))
		print("node is not a goal node")
		local minimum = math.huge
		if goalAdjacent then
			minimum = getTraversalCost(voxelMap, node, goalNode) * heuristic(node.Position, goalNode.Position)
		else
			local faces = getAdjacentFaces(nodeMap, node)
			for index, face in pairs(faces) do
				local cost = computeCost(ship, voxelMap, node, face)
				if cost < minimum then
					print("Found new minimum")
					minimum = cost
					node.ChangeColor(BrickColor.new("Mauve"))
				end
			end
		end
		setLength(ship, node, "LookAhead", minimum)
	end
	--print("remove node:")
	--print("Node Index: ("..node.Index[1]..", "..node.Index[2]..", "..node.Index[3]..")")
	--print("-----")
	queue:Remove(node)
	print("checking node: goal - "..goalLength.." lookahead - "..getLength(ship, node, "LookAhead"))
	if goalLength ~= getLength(ship, node, "LookAhead") then
		print("add node")
		local keys = computeKeys(ship, startNode, node)
		print("keys: ("..keys[1]..", "..keys[2]..")")
		queue:Add(node, keys)
	end
	print("-------------")
end

local function computeShortestPath(ship, nodeMap, voxelMap, queue, startNode, goalNode)
	local topNode, priorityKey = queue:Peek()
	--
	while (comparator(computeKeys(ship, startNode, startNode), priorityKey)) or (getLength(ship, startNode, "LookAhead") ~= getLength(ship, startNode, "Goal")) do
		warn(">>>>>>>>>computing:<<<<<<<<<<<")
		topNode, priorityKey = queue:Pop()--smallest priority
		print("Node Index: ("..topNode.Index[1]..", "..topNode.Index[2]..", "..topNode.Index[3]..")")
		print("Priority Key: ("..priorityKey[1]..", "..priorityKey[2]..")")
		local updatedNodes = {}
		topNode.ChangeColor(BrickColor.new("Really red"))
		local goalLength = getLength(ship, topNode, "Goal")
		local lookAheadLength = getLength(ship, topNode, "LookAhead")
		local goalAdjacent = topNode == goalNode
		
		local newPriorityKey = computeKeys(ship, startNode, topNode)
		if comparator(newPriorityKey, priorityKey) then--if the node's distances have changed then the node is re-added
			queue:Add(topNode, newPriorityKey)
		elseif goalLength > lookAheadLength then
			print("Current goal: "..goalLength)
			print("Current lookahead: "..lookAheadLength)
			setLength(ship, topNode, "Goal", lookAheadLength)
			print("New goal: "..getLength(ship, topNode, "Goal"))
			for index, node in pairs(getNeighboringNodes(nodeMap, topNode)) do
				--wait(0.25)
				updateState(ship, nodeMap, voxelMap, queue, goalNode, startNode, node, goalAdjacent)
				if node.ChangeColor(BrickColor.new("Royal purple")) then
					table.insert(updatedNodes, node)
				end
			end
		else
			print("goal length is less then or equal to lookahead")
			setLength(ship, topNode, "Goal", math.huge)
			for index, node in pairs({topNode, unpack(getNeighboringNodes(nodeMap, topNode))}) do
				--wait(1)
				updateState(ship, nodeMap, voxelMap, queue, goalNode, startNode, node, goalAdjacent)
				if node.ChangeColor(BrickColor.new("Royal purple")) then
					table.insert(updatedNodes, node)
				end
			end
		end
		
		--wait(3)
		for index, node in pairs(updatedNodes) do
			node.ResetColor()
		end
		topNode, priorityKey = queue:Peek()
		queue:Print2(true)
		wait(1)
	end
end

--[[
one can then trace back a shortest path from startNode to any
vertex u by always moving from the current vertex s, starting at u,
to any predecessor s' that minimizes g(s') + c(s', s) until startNode
is reached (ties can be broken arbitrarily)
--]]
local function extractPath()
	local path = {}
	
	return path
end

local function moveShip(ship, path)
	local currentLocation = ship.PrimaryPart.Position 
	for index, position in pairs(path) do
		print("move to point")
		local nextLocation = path[index]
		local movement
		movement = RunService.Heartbeat:Connect(function(step)
			--[[if not pathfinding then
				movement:Disconnect()
				movement = nil
			end]]
			local offsetVector = (nextLocation - currentLocation)
			local distance = offsetVector.magnitude
			if distance > 0 then
				local maxSpeed = SHIP_SPEED < distance and SHIP_SPEED or distance
				local velocity = offsetVector.unit * maxSpeed
				if SHIP_SPEED < distance then
					local shipLocation = currentLocation + (velocity * step)
					ship:SetPrimaryPartCFrame(CFrame.new(shipLocation, nextLocation))-- * CFrame.Angles(0, 90, 0))
					currentLocation = shipLocation
				else
					local shipLocation = ship.PrimaryPart.CFrame + (velocity * step)
					ship:SetPrimaryPartCFrame(shipLocation)
					currentLocation = nextLocation
					movement:Disconnect()
					print("disconnect movement")
					movement = nil
				end
			else
				movement:Disconnect()
				print("disconnect movement")
				movement = nil
			end
		end)
		while movement do
			wait()
		end
		print("stopped moving")
	end
	
end

local function buildFaceVisual(position)
	local visual = Instance.new("Part")
	visual.Anchored = true
	visual.CanCollide = false
	visual.BrickColor = BrickColor.new("Royal purple")
	visual.BottomSurface = Enum.SurfaceType.Smooth
	visual.TopSurface = Enum.SurfaceType.Smooth
	visual.Size = NODE_SIZE
	visual.Position = position
	--visual.Name = "Face"
	--visual.Parent = Workspace
	
	return visual
end

function Pathfinding:FindPath(ship, mapName, currentPosition, targetPosition)

	local queue = priorityQueue.new(comparator)
	local nodeMap = voxelMap.AllNodeMaps[mapName]
	local voxelMap = voxelMap.AllVoxelMaps[mapName]
	if nodeMap and voxelMap then
		local startNode = nodeMap[1][1][1]--getNearestNode(nodeMap, currentPosition)
		--local goalNode = nodeMap[4][4][4]--getNearestNode(nodeMap, targetPosition)
		local goalFace = {nodeMap[3][4][3], nodeMap[4][4][3], nodeMap[3][3][3], nodeMap[4][3][3]}
		local goalNode = goalFace[1]
		
		allLengths[ship] = {}
		local queue = priorityQueue.new(comparator)
		
		startNode.ChangeColor(BrickColor.new("Bright yellow"))
		setLength(ship, startNode, "Goal", math.huge)
		setLength(ship, startNode, "LookAhead", math.huge)
		setLength(ship, goalNode, "Goal", math.huge)
		setLength(ship, goalNode, "LookAhead", 0)
		goalNode.ChangeColor(BrickColor.new("Bright green"))
		queue:Add(goalNode, computeKeys(ship, startNode, goalNode))
		computeShortestPath(ship, nodeMap, voxelMap, queue, startNode, goalNode)
		--[[while pathfinding do
			computeShortestPath(ship, nodeMap, voxelMap, queue, startNode, goalFace)
			--local path = extractPath()
			--moveShip(ship, path)
			--wait for changes in cell traversal cost
			pathfinding = false
		end]]
		
		allLengths[ship] = nil
	else
		warn("Cannot pathfind without nodeMap")
	end
end

--[[
	Visuals:
					 start
			<-----t
	   s0		  s1
		----------  u
		|		 |  |
	    |		 |  |
		|		 |  V
		----------
	   s3         s2
	
		
			t0	  start
		----------
		|		 |
	 u1 |		 | u0
		|		 |
		----------
	goal	t1
	
	to get adjacent faces it can be broken down into getting the four seperaces faces for each of the 6 'big faces'
	
	how will I get the current node's index? (need it to figure out adjacent nodes)
	either:
	-save 'Index' property in nodes
	-pass it through pathfinding (coordinates will probably be retrieved from getNearestNode(position) - in order to get goal and start nodes)
--]]

--[[
	local xOrder = {{0, 1, 0}, {0, 0, 1 * index}, {0, -1, 0}, {0, 0, -1 * index}}

	--must move in both y & z directions - get all 8 nodes
	--must verify that each node exists
	--go in order {y, x * index, -y, -x * index} and return each subsequent node in order
	--once you go through this order, move the first index to the end (go through 4 times)
--]]
--[[
	g(s) -> estimate cost of optimal path from state s to goal
	rhs(s) -> one-step look ahead cost
	
	>get all adjacent faces (24 of them)
	>Faces will be tables containing their 4 vertices)
	>Can use vertices to get each edge
	
	To get coordinates of surface minima(t, u):
	>Minima of an edge: g(se) = y * g(s2) + (1 − y)g(s1)
	>Minima of a face:
		lineOne = ((t1 − t0) * u0 + t0) / (1 − (t1 − t0) * (u1 − u0))
		lineTwo = (u1 − u0) * lineOne + u0
		>minima is the intersection of these two lines
	Using surface minima coords(t, u):
	C -> traversal cost of voxel (that both the start node and the surface share)
	g(surface) = [(g(s1)+(g(s0) − g(s1)) * t] * (1 − u) + [g(s2)+(g(s3) − g(s2)) * t] * u
	rhs(surface) = = C * sqrt(1 + t2 + u2) + [(g(s1)+(g(s0) − g(s1)) * t] * (1 − u) + [g(s2)+(g(s3) − g(s2)) * t] * u
--]]
--check LOS each move, if there is LOS then just move to target
	--otherwise use pathfinding
return Pathfinding
