# Multilevel and Longitudinal Models

Source for an online textbook on multilevel and longitudinal modeling, built with [Quarto](https://quarto.org).

## Building locally

```bash
quarto preview
```

The rendered book goes to `_book/`. Pushing to `main` builds and publishes it to GitHub Pages via `.github/workflows/publish.yml` (set **Settings → Pages → Source** to "GitHub Actions" first).

## Layout

- `_quarto.yml`: book structure and chapter order
- `index.qmd`: preface
- `chapters/`: one `.qmd` per chapter, with figures in `chapters/images/<chapter>/`
- `code/`: R scripts that generate each chapter's figures (run from the repo root)
- `data/`: the datasets used in the book
- `custom.scss`: theme

## License

The text and figures are © 2026 Josh Gilbert and Bill Young under [CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/). The code in `code/` and the code examples in the chapters are under the [MIT License](LICENSE-CODE). The datasets in `data/` keep the terms of their original sources. See [LICENSE](LICENSE) for details.
