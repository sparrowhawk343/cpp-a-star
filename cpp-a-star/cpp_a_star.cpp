#include <algorithm>
#include <map>
#include <unordered_set>
#include <vector>
#include <optional>
#include <queue>

struct Coordinate
{
    int GridX;
    int GridY;

    bool operator!=(const Coordinate& Other) const
    {
        return Tie() != Other.Tie();
    }

    bool operator==(const Coordinate& Other) const
    {
        return Tie() == Other.Tie();
    }

    std::tuple<int, int> Tie() const
    {
        return {GridX, GridY};
    }
};

struct TupleHash
{
    template <typename T, typename U>
    std::size_t operator()(const std::pair<T, U>& Pair) const
    {
        auto hash1 = std::hash<T>{}(Pair.first);
        auto hash2 = std::hash<U>{}(Pair.second);
        return hash1 ^ (hash2 << 1);
    }
};

struct CoordinateHash
{
    std::size_t operator()(const Coordinate& Coord) const
    {
        return TupleHash{}(std::make_pair(Coord.GridX, Coord.GridY));
    }
};

struct CellDefinition
{
    bool IsTraversable;
};

struct CoordinateComparator
{
    bool operator()(const Coordinate A, const Coordinate B) const
    {
        return A.Tie() < B.Tie();
    }
};

struct Cell
{
    Coordinate Coordinates;
    std::optional<Coordinate> ParentCoordinates;
    int GCost, HCost;

    int FCost() const
    {
        return GCost + HCost;
    }

    bool operator==(const Cell& other) const
    {
        return Coordinates == other.Coordinates;
    }
};

struct CellHash
{
    std::size_t operator()(const Cell& cell) const
    {
        const std::size_t Hash = CoordinateHash()(cell.Coordinates);
        return Hash;
    }
};

struct CellComparator
{
    bool operator()(const Cell& A, const Cell& B) const
    {
        return A.FCost() > B.FCost();
    }
};

// because we cannot move diagonally, we'll use the Manhattan Heuristic for our distance check
inline int GetManhattanDistance(const Coordinate& a, const Coordinate& b)
{
    return std::abs(a.GridX - b.GridX) + std::abs(a.GridY - b.GridY);
}

std::vector<Coordinate> GetAccessibleNeighbors(const Coordinate& CellToCheck, const std::vector<std::vector<CellDefinition>>& Grid)
{
    std::vector<Coordinate> AccessibleNeighbors;

    // we get all the neighbors and check if they are traversable
    for (int y = -1; y <= 1; ++y)
    {
        for (int x = -1; x <= 1; ++x)
        {
            // we skip the middle cell (which is CellToCheck) and the diagonal neighbors
            if (x == 0 && y == 0 || x != 0 && y != 0)
            {
                continue;
            }
            const int CheckX = CellToCheck.GridX + x;
            const int CheckY = CellToCheck.GridY + y;

            // out of bounds check for neighbors
            if (CheckY >= 0 && static_cast<std::size_t>(CheckY) < Grid.size() &&
                CheckX >= 0 && static_cast<std::size_t>(CheckX) < Grid[0].size())
            {
                if (Grid[CheckY][CheckX].IsTraversable)
                {
                    AccessibleNeighbors.push_back({CheckX, CheckY});
                }
            }
        }
    }
    return AccessibleNeighbors;
}

std::vector<std::vector<CellDefinition>> InitializeGrid(const std::vector<int>& Map,
                                                        const std::pair<int, int> MapDimensions)
{
    std::vector<std::vector<CellDefinition>> Grid(MapDimensions.second, std::vector<CellDefinition>(MapDimensions.first));
    int CurrentIndex = 0;
    for (int y = 0; y < MapDimensions.second; ++y)
    {
        for (int x = 0; x < MapDimensions.first; ++x)
        {
            // in this case, the implicit conversion of an int of 0 or 1 to bool did not work, so we make it explicit
            const CellDefinition NewDefinition = {Map[CurrentIndex] != 0};
            Grid[y][x] = NewDefinition;
            ++CurrentIndex;
        }
    }
    return Grid;
}

std::vector<int> RetracePath(const Coordinate Start, const Coordinate End, const std::pair<int, int> MapDimensions, const std::map<Coordinate, Cell, CoordinateComparator>& Path)
{
    std::vector<Coordinate> PathAsCoordinates;
    Coordinate Current = End;

    while (Current != Start)
    {
        PathAsCoordinates.push_back(Current);
        Current = Path.at(Current).ParentCoordinates.value();
    }

    std::reverse(PathAsCoordinates.begin(), PathAsCoordinates.end());
    std::vector<int> PathAsIntegers;
    // we know that PathAsIntegers will be as large as PathAsCoordinates, so we reserve space ahead of time
    PathAsIntegers.reserve(PathAsCoordinates.size());

    for (const Coordinate Coord : PathAsCoordinates)
    {
        int Index = Coord.GridY * MapDimensions.first + Coord.GridX;
        PathAsIntegers.push_back(Index);
    }
    return PathAsIntegers;
}

bool FindPath(std::pair<int, int> Start,
              std::pair<int, int> Target,
              const std::vector<int>& Map,
              std::pair<int, int> MapDimensions,
              std::vector<int>& OutPath)
{
    // TODO:
    // bug: when a new shorter path to cell is found, it is not added to PQ, so it's only explored at the score of its more expensive path
    // performance of open set check - guard against fully explored cells to avoid infinite loop - tag popped cells as "do not need exploring"
    // improvement for grid: read up on 2d arrays stored as 1d arrays - std::mdspan from C++23.
    
    const std::vector<std::vector<CellDefinition>> Grid = InitializeGrid(Map, MapDimensions);
    std::priority_queue<Cell, std::vector<Cell>, CellComparator> Open;
    std::unordered_set<Cell, CellHash> OpenSet;
    std::map<Coordinate, Cell, CoordinateComparator> PathSoFar;
    const Coordinate Destination = {Target.first, Target.second};
    const Coordinate StartCoordinate = {Start.first, Start.second};

    bool IsPathFound = false;

    const Cell StartCell{StartCoordinate, std::nullopt, 0, 0};
    Open.emplace(StartCell);
    OpenSet.emplace(StartCell);
    PathSoFar[StartCoordinate] = StartCell;

    while (!Open.empty())
    {
        const Cell Current = Open.top();
        Open.pop();
        OpenSet.erase(Current);

        if (Current.Coordinates == Destination)
        {
            IsPathFound = true;
            OutPath = RetracePath(StartCoordinate, Destination, MapDimensions, PathSoFar);
            break;
        }
        const auto Neighbors = GetAccessibleNeighbors(Current.Coordinates, Grid);
        // go through all accessible neighbors, compute their costs
        for (const Coordinate NeighborCoordinate : Neighbors)
        {
            const int NewMovementCost = PathSoFar.at(Current.Coordinates).GCost + GetManhattanDistance(Current.Coordinates, NeighborCoordinate);
            auto Iterator = PathSoFar.find(NeighborCoordinate);

            if (Iterator == PathSoFar.end() || NewMovementCost < Iterator->second.GCost)
            {
                // only emplace if the PathSoFar map does not yet contain the neighbor cell
                Iterator = PathSoFar.emplace_hint(Iterator, NeighborCoordinate, Cell{NeighborCoordinate, Current.Coordinates, NewMovementCost, GetManhattanDistance(NeighborCoordinate, Destination)});
                // here is where it comes in handy to have this unordered_set that mirrors the PQ - we can easily look up elements in it
                if (std::none_of(OpenSet.begin(), OpenSet.end(), [NeighborCoordinate](const Cell& A) { return A.Coordinates == NeighborCoordinate; }))
                {
                    Open.emplace(Iterator->second);
                    OpenSet.emplace(Iterator->second);
                }
            }
        }
    }

    return IsPathFound;
}

int main(int argc, char* argv[])
{
    // some test maps here:
    // const std::vector<int> Map =
    // {
    //     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    //     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    //     1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
    //     1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1,
    //     1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0,
    //     1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1,
    //     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    //     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    //     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    //     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    //     1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0,
    //     1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1,
    //     1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0,
    //     0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    //     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    //     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    //     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
    // };

     // const std::vector<int> Map = {
     //     1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     //     1, 1, 1, 0, 1, 1, 1, 1, 0, 1,
     //     1, 0, 1, 0, 1, 1, 1, 0, 0, 0,
     //     1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     //     1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     //     1, 1, 1, 1, 1, 1, 1, 1, 0, 1,
     //     1, 0, 0, 0, 1, 1, 0, 1, 0, 0,
     //     1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     //     1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     //     1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
     // };
    const std::vector<int> Map = {1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1};
    // std::vector<int> Map = {0, 0, 1, 0, 1, 1, 1, 0, 1};
    std::vector<int> OutPath;

    std::printf("This is a test print\n");
    // const bool WasPathFound = FindPath({19, 5}, {80, 6}, Map, {100, 100}, OutPath);
     // const bool WasPathFound = FindPath({7, 2}, {9, 0}, Map, {10, 10}, OutPath);
    const bool WasPathFound = FindPath({0, 0}, {1, 2}, Map, {4, 3}, OutPath);
    // const bool WasPathFound = FindPath({2, 0}, {0, 2}, Map, {3, 3}, OutPath);

    for (const int i : OutPath)
    {
        std::printf("i = %d \n", i);
    }

    std::printf("WasPathFound = %d \n", WasPathFound);

    return 0;
}
