-- Cubic Marching Squares algorithm
--
-- Copyright (c) 2012 Reinoud Zandijk (reinoud@13thmonkey.org)
--	All rights reserved.
--
-- Redistribution and use in source and binary forms, with or without
-- modification, are permitted provided that the following conditions
-- are met:
-- 1. Redistributions of source code must retain the above copyright
--    notice, this list of conditions and the following disclaimer.
-- 2. Redistributions in binary form must reproduce the above copyright
--    notice, this list of conditions and the following disclaimer in the
--    documentation and/or other materials provided with the distribution.
--
-- THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND
-- ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
-- IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
-- ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
-- FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
-- DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
-- OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
-- HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
-- LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
-- OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
-- SUCH DAMAGE.
--

module Graphics.Implicit.Export.CubicMarchingSquares (getMesh, gen_octree, calc_corners, enum_edges, interpolate, calc_intersect, gen_subcells, gen_leaf_cells, enum_faces, gen_segments, gen_facelines, trace_obj, trace_objects, triangulation, gen_triangles, sharp_matrix, sharp_inproduct, sharp_aproximation, sharp_center) where

import Graphics.Implicit.Definitions
import Control.Parallel (par, pseq)
import Graphics.Implicit.Export.Util
import qualified Data.List
import qualified Graphics.Implicit.SaneOperators as S

{-
 Cube description:
         3 ________ 2           _____2__     
         /|       /|         / |       /|    
       /  |     /  |      11/  3   10/  |    
   7 /_______ /    |      /__6_|__ /    |1   
    |     |  |6    |     |     |  |     |
    |    0|__|_____|1    |     |__|__0__|    
    |    /   |    /      7   8/   5    /     
    |  /     |  /        |  /     |  /9      
    |/_______|/          |/___4___|/         
   4          5                  

   point 0 is at (x1, y1, z1), point 6 is at (x2, y2, z2)
-}

-- | all calculated info about a point
type Sample = (ℝ3, ℝ3, ℝ)

-- | an Edgecrossing has either no or one crosspoint
type Edgecrossing = [Sample]

-- | a Face : corners, 3D->2D mapping and calculated edge crossings
type Face = ([Sample], ℝ3, ℝ3, ℝ3, [Edgecrossing])

-- | lines on a face have a start, a calculated intermediate point and an end point
type FaceLines = [(Sample, Sample, Sample)]

-- | a Cell has a bounding box, the 8 corners, the 12 edge crossings and for each face the lines
type Cell = (ℝ3, ℝ3, [Sample], [Edgecrossing], [FaceLines])


-- | calculate crossproduct of two vectors
cross :: ℝ3 -> ℝ3 -> ℝ3
cross (v1, v2, v3) (w1, w2, w3) =
	(v2*w3-v3*w2, v3*w1-v1*w3, v1*w2-v2*w1)


-- | calculate a norm for a vertex
normvertex :: ℝ -> Obj3 -> ℝ3 -> (ℝ3, ℝ3)
normvertex res obj p@(x,y,z) = 
	let
		d :: ℝ3 -> ℝ
		d v = ( obj p - obj (p S.+ res S.* v) )
		dx = d (1, 0, 0)
		dy = d (0, 1, 0)
		dz = d (0, 0, 1)
		nonUnitNormal = (dx,dy,dz)
		normal = nonUnitNormal S./ S.norm nonUnitNormal
	in ((x,y,z), normal)

-- | calculate the cube's (end) points, their normals and their values
calc_corners :: ℝ3 -> ℝ3 -> Obj3 -> [Sample]
calc_corners bot top obj =
	let
		eps = 1e-5
		(x1, y1, z1) = bot
		(x2, y2, z2) = top
	in
		[(p, snd $ normvertex eps obj p, obj p) |
			p <- [ 
				(x1, y1, z1), 
				(x2, y1, z1),
				(x2, y2, z1),
				(x1, y2, z1),
				(x1, y1, z2),
				(x2, y1, z2),
				(x2, y2, z2),
				(x1, y2, z2)
			]
		]

-- | create a set of subcells from a specified rectangle and subdevisions
gen_subcells :: ℝ3 -> ℝ3 -> ℝ3 -> [Cell]
gen_subcells bot top divs =
	let
		(x1, y1, z1) = bot
		(x2, y2, z2) = top
		(nx, ny, nz) = divs
		dx = (x2-x1) / nx
		dy = (y2-y1) / ny
		dz = (z2-z1) / nz

		gen :: ℝ3 -> Cell
		gen (xs, ys, zs) = (p0, p1, [], [], [])
			where
				x = x1 + xs*dx
				y = y1 + ys*dy
				z = z1 + zs*dz
				p0 = (x, y, z)
				p1 = (x + dx, y + dy, z + dz)
	in
		[gen (xs, ys, zs) | xs <- [0..nx-1], ys <- [0..ny-1], zs <- [0..nz-1] ]


-- | enumerate all the edges as functions of the points specified
enum_edges :: [Sample] -> [(Sample, Sample)]
enum_edges (p0:p1:p2:p3:p4:p5:p6:p7:[]) =
	[(p0, p1),	--  0
	 (p1, p2),	--  1
	 (p2, p3),	--  2
	 (p3, p0),	--  3
	 (p4, p5),	--  4
	 (p5, p6),	--  5
	 (p6, p7),	--  6
	 (p7, p4),	--  7
	 (p0, p4),	--  8
	 (p1, p5),	--  9
	 (p2, p6),	-- 10
	 (p3, p7)]	-- 11


-- | enumerate all the faces of a cell's samples with corresponding edges per face
enum_faces :: [Sample] -> [Edgecrossing] -> [Face]
enum_faces (p0:p1:p2:p3:p4:p5:p6:p7:[]) (e0:e1:e2:e3:e4:e5:e6:e7:e8:e9:e10:e11:[]) =
	[([p0, p3, p2, p1], x, y, z, [e3, e2, e1, e0]),  -- far
	 ([p4, p5, p6, p7], x, y, z, [e4, e5, e6, e7]),  -- near
	 ([p5, p1, p2, p6], y, z, x, [e9, e1, e10, e5]), -- right
	 ([p0, p4, p7, p3], y, z, x, [e8, e7, e11, e3]), -- left
	 ([p0, p1, p5, p4], x, z, y, [e0, e9, e4, e8]),  -- bottom
	 ([p2, p3, p7, p6], x, z, y, [e2, e11, e6, e10])]-- top
	 	where
			x = (1,0,0)
			y = (0,1,0)
			z = (0,0,1)
	
-- | sort two samples for more deterministic calculation outcomes
sort_sample :: Sample -> Sample -> (Sample, Sample)
sort_sample s1@(p1, n1, v1) s2@(p2, n2, v2) =
	if p1 < p2 then
		(s1, s2)
	else
		(s2, s1)

-- | interpolate between two Samples to find the zero
interpolate ::  ℝ3 -> ℝ -> ℝ3 -> ℝ -> Obj3 -> (ℝ, ℝ3)
interpolate left@(x1, y1, z1) left_obj right@(x2, y2, z2) right_obj obj =
	let
		eps = 1e-2
		mid = (x1 + hdx, y1 + hdy, z1 + hdz);
		mid_obj = obj mid
		hdx = (x2 - x1)/2.0
		hdy = (y2 - y1)/2.0
		hdz = (z2 - z1)/2.0
	in
		if (hdx < eps) && (hdy < eps) && (hdz < eps) then
			(0, mid)
		else
			if (left_obj * mid_obj) < 0 then
				interpolate left left_obj mid mid_obj obj
			else
				interpolate mid mid_obj right right_obj obj

-- | calculate the intersection if present of the edge with the object
calc_intersect_edge :: Obj3 -> (Sample, Sample) -> [Sample]
calc_intersect_edge obj (s1, s2) =
	let
		eps = 1e-5
		((p1, _, v1), (p2, _, v2)) = sort_sample s1 s2
		(guess, hit) = interpolate p1 v1 p2 v2 obj
		hitv = obj hit
		res = (hit, snd $ normvertex eps obj hit, hitv)
	in
		if (obj p1)*(obj p2) > 0 then
			[]
		else
			if (guess >= 0) || (guess <= 1) then 
				[res]
			else
				[]

-- | returns if there are corners that are too sharp
needs_split :: [Sample] -> [Edgecrossing] -> Bool
needs_split corners edges = not. null $ filter (<0.7) (map calc_angles $ enum_edges corners)
	where
		calc_angles ((_, n1, _), (_, n2, _)) = n1 S.⋅ n2

-- | create a signed adaptive octree
gen_octree :: ℝ3 -> ℝ3 -> ℝ -> Obj3 -> [Cell]
gen_octree bot top level obj =
	let
		(x1, y1, z1) = bot
		(x2, y2, z2) = top
		subcells = gen_subcells bot top (2,2,2)
		corners = calc_corners bot top obj
		edges = map (calc_intersect_edge obj) (enum_edges corners)
	in
		if (null $ concat edges) then
			[]
		else
		if (level > 0) {- && needs_split corners edges -} then
			-- create our subcells
			concat $ [gen_octree sbot stop (level-1) obj |
				(sbot, stop, _, _, _) <- subcells]
		else
			-- generate our output
			[(bot, top, corners, edges, [])]

-- | generate all the cells in the specified box
gen_leaf_cells :: ℝ3 -> ℝ3 -> ℝ -> Obj3 -> [Cell]
gen_leaf_cells p@(x1, y1, z1) q@(x2, y2, z2) res obj = 
	let
		p'@(x1', y1', z1') = (x1-1.0, y1-1.0, z1-1.0)
		q'@(x2', y2', z2') = (x2+1.0, y2+1.0, z2+1.0)

		-- How many steps will we take on each axis?
		nx = fromIntegral $ ceiling $ (x2' - x1') / res
		ny = fromIntegral $ ceiling $ (y2' - y1') / res
		nz = fromIntegral $ ceiling $ (z2' - z1') / res

		-- fixed max depth for now
		cells = [gen_octree b p 3 obj | (b, p, _, _, _) <- gen_subcells p' q' (nx, ny, nz)]
	in
		[f | f@(b, t, c, e, l) <- concat cells, not . null $ concat e]



-- | calculate intersection between two samples using their normals
calc_intersect :: ℝ3 -> ℝ3 -> Sample -> Sample -> ℝ3 -> ℝ3 -> ℝ3 -> ℝ3
calc_intersect bot top s0 s1 eA eB eC =
	let
		((p0, n0, v0), (p1, n1, v1)) = sort_sample s0 s1
		eps = 1e-5
		-- get the tangents
		n0x = eA S.⋅ n0
		n0y = eB S.⋅ n0
		n1x = eA S.⋅ n1
		n1y = eB S.⋅ n1

		v0x = -n0y
		v0y =  n0x
		v1x = -n1y
		v1y =  n1x

		-- C axis is equal for p0 and p1
		pz = eC S.⋅ p0

		-- generic 2 line crosser
		x1 = eA S.⋅ p0
		y1 = eB S.⋅ p0
		x2 = x1 + v0x
		y2 = y1 + v0y
		x3 = eA S.⋅ p1
		y3 = eB S.⋅ p1
		x4 = x3 + v1x
		y4 = y3 + v1y

		xmin = eA S.⋅ bot
		xmax = eA S.⋅ top
		ymin = eB S.⋅ bot
		ymax = eB S.⋅ top

		denom  = (y4-y3) * (x2-x1) - (x4-x3) * (y2-y1)
		numera = (x4-x3) * (y1-y3) - (y4-y3) * (x1-x3)
		numerb = (x2-x1) * (y1-y3) - (y2-y1) * (x1-x3)

		mua = numera / denom
--		mub = numerb / denom
		xs = x1 + mua * (x2 - x1)
		ys = y1 + mua * (y2 - y1)

		mid = (((x1 + x3)/2.0) S.* eA) S.+ (((y1 + y3)/2.0) S.* eB) S.+ (pz S.* eC)
		parallel  = (abs denom < eps) || ((abs numera) < eps && (abs numerb) < eps)
		too_sharp = xs <= xmin || xs >= xmax || ys <= ymin || ys >= ymax 
	in
		if too_sharp || parallel then
			-- when parallel, or has too sharp angles, just take mid point
			mid
		else
			-- return calculated intersection
			xs S.* eA S.+ ys S.* eB S.+ pz S.* eC


-- | generate the slice lines for a face on a cube
gen_facelines :: ℝ3 -> ℝ3 -> Face -> FaceLines
gen_facelines bot top (samples, eA, eB, eC, edgecrossings) =
	let
		s0:s1:s2:s3:[] = samples
		ec0:ec1:ec2:ec3:[] = edgecrossings
		(p0,_,p0v) = s0
		(p1,_,p1v) = s1
		(p2,_,p2v) = s2
		(p3,_,p3v) = s3
		e0 = head ec0
		e1 = head ec1
		e2 = head ec2
		e3 = head ec3

		hat :: Sample -> Sample -> FaceLines
		hat a b = [(a, ab, b)]
			where
				intersection = calc_intersect bot top a b eA eB eC
				ab = (intersection, (0,0,0), 0)
	in
		case (p0v < 0, p1v < 0, p2v < 0, p3v < 0) of
		(False, False, False, False) ->	-- not possible
			[]
		(True,  True,  True,  True)  -> -- not possible?
			[]
		-- singletons
		(True,  False, False, False) -> -- p0v < 0, rest positive
			hat e0 e3 
		(False, True,  False, False) -> -- p1v < 0, rest positive
			hat e0 e1 
		(False, False, True,  False) -> -- p2v < 0, rest positive
			hat e1 e2 
		(False, False, False, True)  -> -- p3v < 0, rest positive
			hat e3 e2 
		-- singletons
		(False, True,  True,  True)  -> -- p0v > 0, rest negative
			hat e3 e0 
		(True,  False, True,  True)  -> -- p1v > 0, rest negative
			hat e1 e0 
		(True,  True,  False, True)  -> -- p2v > 0, rest negative
			hat e2 e1 
		(True,  True,  True,  False) -> -- p3v > 0, rest negative
			hat e2 e3 
		-- vertical and horizontal deviders
		(True,  True,  False, False) -> -- p0v and p1v negative
			hat e1 e3 
		(False, False, True, True)   -> -- p2v and p3v negative
			hat e3 e1
		(True,  False, False, True)  -> -- p0v and p3v negative
			hat e0 e2 
		(False, True,  True,  False) -> -- p1v and p2v negative
			hat e2 e0 
		-- ambiguos cased (p0v, p2v) or (p1v, p3v)
		(True,  False, True, False)  -> -- p0v and p2v negative
			error $ "Ambiguous case tftf"
--			[]	-- XXX TODO resolve me
		(False, True,  False, True)  -> -- p1v and p3v negative
			error $ "Ambiguous case ftft"
--			[]	-- XXX TODO resolve me
			

-- | create object segments for the cell by calculating the face lines
gen_segments :: Cell -> Cell
gen_segments cell@(bot, top, corners, edges, _) =
	(bot, top, corners, edges, map (gen_facelines bot top) $ enum_faces corners edges)


-- | trace an object in a cell by picking the lines (from all the faces) that form a cycle.
trace_obj :: (Sample, Sample, Sample) -> (FaceLines, FaceLines) -> (FaceLines, FaceLines)
trace_obj start (object, []) = (object, [])
trace_obj start (object, (l:rest)) =
	let
		isend (_,_,p1) (p2,_,_) = (p1 == p2)
		prepend (a,_,b) (c,_,d) = (d == a)
		swap (a,i,b) = (b,i,a)
		sw_l = swap l
	in
		if prepend (head object) l then
			if (isend start l) then
				(l:object, rest)
			else
				trace_obj start (l:object, rest)
		else
		if prepend (head object) sw_l then
			if (isend start sw_l) then
				(sw_l:object, rest)
			else
				trace_obj start (sw_l:object, rest)
		else
			trace_obj start (object, rest ++ [l])


-- | trace all objects (cycles) from the lines in a cell
trace_objects :: [FaceLines] -> FaceLines -> [FaceLines]
trace_objects done [] = done
trace_objects done (start:rest) =
	let
		(object, left) = trace_obj start ([start], rest)
	in
		trace_objects (object:done) left
--		object:done	-- IFF you only want to see the 1st object


-- | get the list of normals in a matrix
sharp_matrix :: FaceLines -> [ℝ3]
sharp_matrix lines = [n | ((_,n,_),_,_) <- lines]


-- | get the inproduct of the points and their normals
sharp_inproduct :: FaceLines -> [ℝ]
sharp_inproduct lines = [p S.⋅ n | ((p,n,_),_,_) <- lines]


-- | calc aproximation
sharp_aproximation :: [ℝ3] -> ℝ3 -> [ℝ]
sharp_aproximation mat vec = map (vec S.⋅) mat


-- | calculate center of lines, ignoring the calculated intersection for now
sharp_center :: FaceLines -> ℝ3
sharp_center lines =
	let
		sum :: FaceLines -> ℝ3
		sum [] = (0,0,0)
		sum (((s1,_,_), (s2,_,_), (s3,_,_)):rest) =
			s1 S.+ s3 S.+ (sum rest)
--			s1 S.+ s2 S.+ s3 S.+ (sum rest)
		(sx, sy, sz) = sum lines
		l = fromIntegral $ 2 * length lines
	in
		(sx/l, sy/l, sz/l)


-- | find best fit of an extremity (cap) on the 3d object defined by
-- the facelines using their normals using least squares method
sharp_find :: FaceLines -> ℝ3
sharp_find lines =
	let
		e = 1e-2
		mat = sharp_matrix lines
		rhs = sharp_inproduct lines
		center = sharp_center lines

		sharp_calc point = (value, point)
			where
				sqd (a,b) = (a - b)**2
				value = sum $ map sqd $ zip (sharp_aproximation mat point) rhs

		sharp_find_iter :: (ℝ, ℝ3) -> ℝ3
		sharp_find_iter (curval, (x, y, z)) =
			let
				env = [sharp_calc (x', y', z') |
					x' <- [x-e, x, x+e],
					y' <- [y-e, y, y+e],
					z' <- [z-e, z, z+e]]
				minpoint = head $ Data.List.sort env
			in
				if (fst minpoint) < curval then
					sharp_find_iter minpoint
				else
					(x, y, z)

		calc_angles ((_, n1, _), (_,_,_), (_, n2, _)) = n1 S.⋅ n2
	in
		-- first detect sharp features before racing outside the box,
		-- this prevents runaways of the cap
		if not . null $ filter (<0.999) (map calc_angles lines) then
			sharp_find_iter $ sharp_calc center
		else
			center

-- | triangulate the lines
triangulation :: FaceLines -> TriangleMesh
triangulation lines =
	let
		sharp = sharp_find lines
		triangles f i e = [(f, sharp, i), (i, sharp, e)]
	in
		concat [triangles f i e | ((f,_,_),(i,_,_),(e,_,_)) <- lines]


-- | explicit spacing of the cubes for visualisation / debug purposes
spacify :: ℝ3 -> TriangleMesh -> TriangleMesh
spacify (x0, y0, z0) mesh = map movetriangle mesh
	where
		movetriangle (a, b, c) = (a S.+ delta, b S.+ delta, c S.+ delta)
		delta = (0.15 * x0, 0.15 * y0, 0.15 * z0)


-- | create triangles for a given give
gen_triangles :: Cell -> TriangleMesh
gen_triangles cell@(bot, top, corners, edges, facelines) =
	let
		objects = trace_objects [] $ concat facelines
		-- TODO de_ambiguatify two caps, rare though
		process obj = {- spacify bot $ -} triangulation obj
	in
		concat $ map process objects


-- | getMesh gets a triangle mesh describe the boundary of your 3D
--  object. 
getMesh :: ℝ3 -> ℝ3 -> ℝ -> Obj3 -> TriangleMesh
getMesh p q res obj = 
	let
		-- step 1, generate leaf cells
		step1 = gen_leaf_cells p q res obj

		-- step 2 generate segments
		step2 = map (gen_segments) step1

		-- step 3 extract surface
		triangles = map (gen_triangles) step2

		-- Remove degenerated triangles
		sane :: Triangle -> Bool
		sane (a, b, c) = not (a == b || a == c || b == c)
	in
		concat $ triangles
		-- still needed?
--		[triangle | triangle <- concat $ triangles, sane triangle]

