--[[
tabset-beamer.lua — rendu Beamer des slides Quarto pensées pour revealjs

À utiliser dans le YAML, pour le format beamer uniquement :

  format:
    beamer:
      include-in-header: mpe-beamer-theme.tex
      filters:
        - at: pre-ast
          path: tabset-beamer.lua

Ce que fait le filtre (uniquement en sortie beamer, sans effet ailleurs) :
  * ::: panel-tabset  -> une frame par onglet. Chaque frame reprend le titre
    de la slide et affiche une barre d'onglets avec l'onglet actif en
    surbrillance (\qtab / \qtabactive, définis dans le préambule).
    Si la slide ne contient rien d'autre que le tabset, le premier onglet
    reste sur la frame d'origine. Pour forcer ce comportement même quand
    du contenu précède le tabset :  ::: {.panel-tabset same-frame="true"}
  * ::: callout-*     -> environnement mpecallout (tcolorbox aux couleurs
    du CSS : fond gris clair, liseré gauche coloré, titre coloré).
  * .column isolé    -> enveloppé automatiquement dans un ::: columns
    (sinon erreur LaTeX).
  * \begin{align} ... -> align* (pas de numéros, comme MathJax).
  * [texte]{.rouge}  -> \textcolor{rouge}{texte}, pour toutes les classes
    de couleur du CSS (jaune, jauneb, orange, rouge, vert, bleu, bleuf,
    clair, alea, fixe, care, highlight).

Contrainte : un tabset doit être au niveau principal de la slide (pas à
l'intérieur d'une colonne ou d'un autre div), puisqu'il crée des frames.
--]]

local function is_beamer()
  if quarto and quarto.doc and quarto.doc.is_format then
    return quarto.doc.is_format("beamer")
  end
  return FORMAT == "beamer"
end

local function slide_level()
  local sl = PANDOC_WRITER_OPTIONS and PANDOC_WRITER_OPTIONS.slide_level
  return tonumber(sl) or 2
end

-- Inlines -> chaîne LaTeX
local function to_latex(inlines)
  if not inlines or #inlines == 0 then return "" end
  local s = pandoc.write(pandoc.Pandoc({ pandoc.Plain(inlines) }), "latex")
  return (s:gsub("%s+$", ""))
end

------------------------------------------------------------------------
-- Couleurs (classes CSS -> \textcolor)
------------------------------------------------------------------------
local color_classes = {
  "jaune", "jauneb", "orange", "rouge", "vert", "bleu", "bleuf",
  "clair", "alea", "fixe", "care", "highlight",
}

local function Span(el)
  for _, c in ipairs(color_classes) do
    if el.classes:includes(c) then
      local res = pandoc.Inlines({ pandoc.RawInline("latex", "\\textcolor{" .. c .. "}{") })
      if c == "care" or c == "highlight" then
        res = pandoc.Inlines({ pandoc.RawInline("latex", "\\textcolor{" .. c .. "}{\\bfseries ") })
      end
      res:extend(el.content)
      res:insert(pandoc.RawInline("latex", "}"))
      return res
    end
  end
end

------------------------------------------------------------------------
-- Callouts -> mpecallout
------------------------------------------------------------------------
local callout_colors = {
  ["callout-note"] = "rouge",
  ["callout-important"] = "rouge",
  ["callout-tip"] = "vert",
  ["callout-warning"] = "orange",
  ["callout-caution"] = "jauneb",
}

local function callout(div)
  -- ::: hidden (macros MathJax pour revealjs) : supprimé en beamer.
  -- Les macros sont déjà fournies par mpe-beamer-theme.tex, et ce bloc
  -- finirait dans une frame où les #1 des \newcommand font planter LaTeX
  -- ("Illegal parameter number in definition of \iterate").
  if div.classes:includes("hidden") then return {} end
  local col
  for cls, c in pairs(callout_colors) do
    if div.classes:includes(cls) then col = c end
  end
  if not col then return nil end

  local content = div.content
  local title = pandoc.Inlines({})
  if div.attributes["title"] then
    local parsed = pandoc.read(div.attributes["title"], "markdown").blocks
    if #parsed > 0 then title = parsed[1].content end
  end
  if #content > 0 and content[1].t == "Header" then
    title = content[1].content
    content = pandoc.Blocks({ table.unpack(content, 2) })
  end

  local out = pandoc.Blocks({
    pandoc.RawBlock("latex",
      "\\begin{mpecallout}[" .. col .. "]{" .. to_latex(title) .. "}"),
  })
  out:extend(content)
  out:insert(pandoc.RawBlock("latex", "\\end{mpecallout}"))
  return out
end

------------------------------------------------------------------------
-- Colonnes isolées -> envelopper dans ::: columns
------------------------------------------------------------------------
local fix_columns

fix_columns = function(blocks)
  local out = pandoc.Blocks({})
  local run = nil
  local function flush()
    if run then
      out:insert(pandoc.Div(run, pandoc.Attr("", { "columns" })))
      run = nil
    end
  end
  for _, b in ipairs(blocks) do
    if b.t == "Div" and b.classes:includes("column") then
      b.content = fix_columns(b.content)
      run = run or pandoc.Blocks({})
      run:insert(b)
    else
      flush()
      if b.t == "Div" then
        if b.classes:includes("columns") then
          -- largeurs sommant à 100 % : réduites à 96 % (place pour l'espacement beamer)
          local total = 0
          for _, c in ipairs(b.content) do
            if c.t == "Div" and c.attributes.width then
              total = total + (tonumber((c.attributes.width:gsub("%%", ""))) or 0)
            end
          end
          for _, c in ipairs(b.content) do
            if c.t == "Div" then
              c.content = fix_columns(c.content)
              local w = c.attributes.width and tonumber((c.attributes.width:gsub("%%", "")))
              if w and total > 96 then
                c.attributes.width = string.format("%.1f%%", w * 96 / total)
              end
            end
          end
        else
          b.content = fix_columns(b.content)
        end
      end
      out:insert(b)
    end
  end
  flush()
  return out
end

------------------------------------------------------------------------
-- Tabsets -> une frame par onglet
------------------------------------------------------------------------
local function parse_tabs(div)
  local level
  for _, b in ipairs(div.content) do
    if b.t == "Header" then level = b.level; break end
  end
  local tabs = {}
  if not level then
    return { { title = pandoc.Inlines({}), content = div.content } }
  end
  for _, b in ipairs(div.content) do
    if b.t == "Header" and b.level == level then
      table.insert(tabs, { title = b.content, content = pandoc.Blocks({}) })
    elseif #tabs > 0 then
      tabs[#tabs].content:insert(b)
    end
  end
  return tabs
end

local function tabbar(tabs, active)
  local parts = {}
  for k, t in ipairs(tabs) do
    local cmd = (k == active) and "\\qtabactive" or "\\qtab"
    table.insert(parts, cmd .. "{" .. to_latex(t.title) .. "}")
  end
  return pandoc.RawBlock("latex",
    "\\begin{qtabbar}" .. table.concat(parts, "\\qtabsep") .. "\\end{qtabbar}")
end

local function split_tabsets(blocks)
  local sl = slide_level()
  local out = pandoc.Blocks({})
  local cur_header = nil
  local since_header = 0
  local counter = 0

  for _, b in ipairs(blocks) do
    if b.t == "Header" and b.level <= sl then
      cur_header = (b.level == sl) and b or nil
      since_header = 0
      out:insert(b)
    elseif b.t == "Div" and b.classes:includes("panel-tabset") then
      counter = counter + 1
      local tabs = parse_tabs(b)
      local same = (b.attributes["same-frame"] == "true") or since_header == 0
      local base_id = (cur_header and cur_header.identifier ~= "")
        and cur_header.identifier or ("tabset-" .. counter)
      for k, tab in ipairs(tabs) do
        if k > 1 or not same then
          local h
          if cur_header then
            h = cur_header:clone()
          else
            h = pandoc.Header(sl, pandoc.Inlines({}))
          end
          h.identifier = base_id .. "-tab" .. k
          out:insert(h)
        end
        out:insert(tabbar(tabs, k))
        out:extend(tab.content)
      end
      since_header = 1
    else
      out:insert(b)
      since_header = since_header + 1
    end
  end
  return out
end

------------------------------------------------------------------------
-- \begin{align} -> \begin{align*} (MathJax ne numérote pas, LaTeX si)
------------------------------------------------------------------------
-- flextable : définition de \ascline (avec #1, #2, #3) émise dans chaque
-- tableau ; dans une frame beamer elle provoque "Illegal parameter number
-- in definition of \iterate". On la retire (elle est dans le thème).
local ascline_def = "\\providecommand{\\ascline}[3]{\\noalign{\\global\\arrayrulewidth #1}\\arrayrulecolor[HTML]{#2}\\cline{#3}}"

local function remove_plain(text, s)
  local i, j = text:find(s, 1, true)
  while i do
    text = text:sub(1, i - 1) .. text:sub(j + 1)
    i, j = text:find(s, 1, true)
  end
  return text
end

local function unnumber(el)
  if el.format == "latex" or el.format == "tex" then
    el.text = remove_plain(el.text, ascline_def)
    el.text = el.text
      :gsub("\\begin{(%a+)}", function(env)
        if env == "align" or env == "equation" or env == "gather" or env == "multline" then
          return "\\begin{" .. env .. "*}"
        end
      end)
      :gsub("\\end{(%a+)}", function(env)
        if env == "align" or env == "equation" or env == "gather" or env == "multline" then
          return "\\end{" .. env .. "*}"
        end
      end)
      :gsub("\\\\%s*(\\end{align%*})", "%1")   -- \\ final inutile
    return el
  end
end

------------------------------------------------------------------------
function Pandoc(doc)
  if not is_beamer() then return nil end
  doc = doc:walk({ Span = Span, Div = callout, RawBlock = unnumber, RawInline = unnumber })
  doc.blocks = fix_columns(doc.blocks)
  doc.blocks = split_tabsets(doc.blocks)
  return doc
end
