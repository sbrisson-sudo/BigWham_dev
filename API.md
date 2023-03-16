## Mermaid charts
https://mermaid.js.org/syntax/classDiagram.html


```mermaid
---
title: BigWham API
---
classDiagram
    BoundaryElement <|-- Segment
    BoundaryElement <|-- Polygon
    Polygon <|-- Triangle
    Polygon <|-- Rectangle

    
    class BoundaryElement{
        +getCollocationPoints()
    }
    class Segment{
    }
    class Triangle{
    }
    class Polygon{
    }
```
