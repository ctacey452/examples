class VertexPair
{
public:
    VertexPair (const TopoDS_Vertex& theFirstVertex, const TopoDS_Vertex& theSecondVertex) :
        myFirstVertex (theFirstVertex), mySecondVertex (theSecondVertex) {}

    bool operator== (const VertexPair& thePair) const
    {
        bool aFirstComparison = thePair.FirstVertex().IsSame (myFirstVertex) && 
            thePair.SecondVertex().IsSame (mySecondVertex);
        bool aSecondComparison = thePair.FirstVertex().IsSame (mySecondVertex) && 
            thePair.SecondVertex().IsSame (myFirstVertex);

        return (aFirstComparison || aSecondComparison);
    }

private:
    TopoDS_Vertex myFirstVertex;
    TopoDS_Vertex mySecondVertex;
};

class HashFunction
{
public:
    //XOR-hash is used intentionally to avoid difficulties with vertex-pair ordering
    size_t operator() (const VertexPair& thePair) const
    {
        return Hash <TopoDS_Vertex>() (thePair.FirstVertex()) ^ Hash <TopoDS_Vertex>() (thePair.SecondVertex());
    }
};

bool CompareEdgesMidpoints (const TopoDS_Edge& theEdge1, const TopoDS_Edge& theEdge2)
{
    auto aMiddlePoint = [] (const TopoDS_Edge& theEdge) {
        double aFirstParameter, aLastParameter;
        TopLoc_Location aLocation;
        auto aCurve = BRep_Tool::Curve (theEdge, aLocation, aFirstParameter, aLastParameter);

        return aCurve->Value ((aFirstParameter + aLastParameter) / 2);
    };

    if (theEdge1.IsNull() || theEdge2.IsNull()) {
        return false;
    }

    return !!aMiddlePoint (theEdge1).IsEqual (aMiddlePoint (theEdge2), 1e-7);
}

bool WireInFace (const TopoDS_Wire& theWire, const TopoDS_Face& theFace)
{
    auto anOuterWire = BRepTools::OuterWire (theFace);
    for (TopExp_Explorer anEdgeExp1 (theWire, TopAbs_EDGE); anEdgeExp1.More(); anEdgeExp1.Next()) {
        if (!EdgeInWire (TopoDS::Edge (anEdgeExp1.Current()), anOuterWire)) {
            return false;
        }
    }
    return true;
}

bool EdgeInWire (const TopoDS_Edge& theEdge, const TopoDS_Wire& theWire)
{
    for (TopExp_Explorer anEdgeExp (theWire, TopAbs_EDGE); anEdgeExp.More(); anEdgeExp.Next()) {
        if (anEdgeExp.Current().IsSame (theEdge)) {
            return true;
        }
    }
    return false;
}

bool ExistCommonEdgesAroundTheHole (
    const TopoDS_Face& theFace1, const TopoDS_Face& theFace2,
    const TopoDS_Wire& theWire1, const TopoDS_Wire& theWire2)
{
    std::array <TopoDS_Face, 2> aFaces = {{ theFace1, theFace2 }};
    std::array <TopoDS_Wire, 2> aWires = {{ theWire1, theWire2 }};

    std::array <TopoDS_Edge, 4> anEdgeVector;
    anEdgeVector.fill (TopoDS_Edge());

    TopoDS_Vertex aVertex1, aVertex2;
    TopExp::Vertices (theWire1, aVertex1, aVertex2);

    for (auto i = 0; i < 2; ++i) {
        auto anOuterWire = BRepTools::OuterWire (aFaces[i]);

        for (TopExp_Explorer anEdgeExp (anOuterWire, TopAbs_EDGE); anEdgeExp.More(); anEdgeExp.Next()) {
            auto anEdge = TopoDS::Edge (anEdgeExp.Current());

            if (EdgeInWire (anEdge, aWires[i])) {
                continue;
            }

            TopoDS_Vertex aFirstPoint, aSecondPoint;
            TopExp::Vertices (anEdge, aFirstPoint, aSecondPoint);

            if (aFirstPoint.IsSame (aVertex1) || aSecondPoint.IsSame (aVertex1)) {
                anEdgeVector[2 * i] = anEdge;
            } else if (aFirstPoint.IsSame (aVertex2) || aSecondPoint.IsSame (aVertex2)) {
                anEdgeVector[2 * i + 1] = anEdge;
            }
        }
        if (anEdgeVector[2 * i].IsNull() || anEdgeVector[2 * i + 1].IsNull()) {
            return false;
        }
    }

    if (anEdgeVector[0].IsSame (anEdgeVector[2]) && anEdgeVector[1].IsSame (anEdgeVector[3])) {
        return true;
    }
    return false;
}

//In a map for each vector of wires for each face try to found wire in another face with the same
//verticies. Its a good way to detect multiface hole. After a hole was found a new edge for intersection 
//line was created and new face for each part of hole was generated.
//All of this are saved in a HoleStructure
std::vector<HoleStructure> FillHoleStructure (
    std::unordered_map <TopoDS_Face, std::vector<TopoDS_Wire>>& theFWMap)
{
    std::vector<HoleStructure> aHoleStructure;
    std::unordered_map <VertexPair,
        std::vector <std::pair <TopoDS_Face, TopoDS_Wire>, HashFunction> aVertexPairFaceMap;

    for (auto anIterMap = theFWMap.begin(); anIterMap != theFWMap.end(); ++anIterMap) {
        auto& aWireVector = anIterMap->second;
        for (auto anIterWires = aWireVector.begin(); anIterWires != aWireVector.end(); ++anIterWires) {
            TopoDS_Vertex aFirstPoint, aSecondPoint;
            TopExp::Vertices (*anIterWires, aFirstPoint, aSecondPoint);

            auto aPair = VertexPair (aFirstPoint, aSecondPoint);
            aVertexPairFaceMap [aPair].push_back (std::make_pair (anIterMap->first, *anIterWires));
        }
    }
    for (auto anIterMap = aVertexPairFaceMap.begin(); anIterMap != aVertexPairFaceMap.end(); ++anIterMap) {
        if (anIterMap->second.size() == 1) {
            continue;
        }

        for (auto anIterFirst = anIterMap->second.begin(); anIterFirst != anIterMap->second.end(); ++anIterFirst) {
            auto aFirstFace = anIterFirst->first;
            auto aFirstWire = anIterFirst->second;

            TopoDS_Vertex aFirstPointFirst, aSecondPointFirst;
            TopExp::Vertices (aFirstWire, aFirstPointFirst, aSecondPointFirst);

            if (aFirstPointFirst.IsSame (aSecondPointFirst)) {
                continue;
            }

            for (auto anIterSecond = std::next (anIterFirst); anIterSecond != anIterMap->second.end(); ++anIterSecond) {
                auto aSecondFace = anIterSecond->first;
                auto aSecondWire = anIterSecond->second;

                if (WireInFace (aFirstWire, aSecondFace)) {
                    continue;
                }

                if (WireInFace (aSecondWire, aFirstFace)) {
                    continue;
                }

                if (!ExistCommonEdgesAroundTheHole (aFirstFace, aSecondFace, aFirstWire, aSecondWire)) {
                    continue;
                }

                TopoDS_Edge anEdge = CreateCommonEdge (aFirstFace, aSecondFace, aFirstPointFirst, aSecondPointFirst);

                if (anEdge.IsNull()) {
                    anEdge = EdgeFromClosingParamSegment (aFirstFace, aFirstWire);
                    TopoDS_Edge aTempEdge = EdgeFromClosingParamSegment (aSecondFace, aSecondWire);

                    if (!CompareEdgesMidpoints (anEdge, aTempEdge)) {
                        anEdge.Nullify();
                    }
                }

                if (!anEdge.IsNull()) {
                    auto AddToHoleStructure = [&] (const TopoDS_Wire& theWire,
                        const TopoDS_Wire& thePartWire,
                        const TopoDS_Face& theBaseFace,
                        const TopoDS_Face& theFirstFace,
                        const TopoDS_Face& theSecondFace) {

                        std::vector<TopoDS_Face> aNewFaces;
                        aNewFaces.push_back (theFirstFace);
                        aNewFaces.push_back (theSecondFace);
                        aHoleStructure.push_back (HoleStructure (theWire, thePartWire, theBaseFace, aNewFaces));
                    };
                    auto aWire = MakeSingleWire (aFirstWire, aSecondWire);
                    auto aNewFirstFace = CreateNewFace (aFirstFace, aFirstWire, anEdge);
                    auto aNewSecondFace = CreateNewFace (aSecondFace, aSecondWire, anEdge);

                    AddToHoleStructure (aWire, aFirstWire, aFirstFace, aNewFirstFace, aNewSecondFace);
                    AddToHoleStructure (aWire, aSecondWire, aSecondFace, aNewSecondFace, aNewFirstFace);
                }
            }
        }
    }

    return aHoleStructure;
}