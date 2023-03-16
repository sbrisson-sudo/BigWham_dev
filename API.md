## Mermaid charts
https://mermaid.js.org/syntax/classDiagram.html


```mermaid
---
title: BigWham API
---
classDiagram
    class BoundaryElement~int dim, int nvert, int p~{
        +getCollocationPoints()
    }
    class Segment{
    }
    class Triangle{
    }
    class Polygon{
    }
    BoundaryElement <|-- Segment
    BoundaryElement <|-- Polygon
    Polygon <|-- Triangle
    Polygon <|-- Rectangle

```
