## Documentation Generation
- We are using [mkdocs](https://www.mkdocs.org/) to generate the documentation. We use specifically the [material theme](https://squidfunk.github.io/mkdocs-material/). You can install the dependencies with the following command:
```bash
pip install mkdocs mkdocs-material
```
- You need to edit `mkdocs.yml` to add new pages.
- To make changes in documentation, edit the markdown files in the `docs` folder. To generate the documentation, run the following commands:
```bash
# to serve locally (required for testing)
mkdocs serve
# to publish using github pages
mkdocs gh-deploy
```

### Testing
I am testing [katex](https://squidfunk.github.io/mkdocs-material/reference/math/?h=tex#katex) here $\mathbf{G}$


$$
\mathbf{t}(\bf{x}) - \mathbf{t}_0(\bf{x}) = \int_S \mathbf{H}(\bf{x} - \bf{y}) \cdot \mathbf{u}(\bf{y}) \, d\bf{y}
$$