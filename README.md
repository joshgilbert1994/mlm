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
- `custom.scss`: theme
