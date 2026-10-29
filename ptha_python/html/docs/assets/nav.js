// Shared navigation of the ptha_python documentation: the page list on
// the left, "On this page" on the right (built from the page's own h2), the
// previous/next links, copy buttons on code blocks, and tab groups.
(function () {
  var GROUPS = [
    ["Start", [
      ["index.html", "Home"],
      ["quickstart.html", "Quick start"],
      ["concepts.html", "PTHA in plain words"],
      ["pipeline.html", "The nine steps"]]],
    ["Geometry: from SLAB to mesh", [
      ["step1.html", "1. Contours from SLAB"],
      ["clean_ends.html", "1. Clean ends"],
      ["step2.html", "2. The mesh"],
      ["discretizer.html", "2. optimal, lm and mid"],
      ["cell_size.html", "Cell size (strasser)"],
      ["mesh_file.html", "Your own mesh (mesh_file)"]]],
    ["Rates", [
      ["step3.html", "3. Plate convergence"],
      ["convergence_sources.html", "3. Bird or Bird + Griffin, per zone"],
      ["seismicity.html", "4-5. Earthquake catalogue"],
      ["input_json.html", "6. The input file"],
      ["engine.html", "7. Logic tree and rates"],
      ["rupture_size.html", "Rupture size"],
      ["segmentation.html", "Segmented vs unsegmented"],
      ["hs_vaus.html", "7b. HS and VAUS slip"],
      ["variable_mu.html#how", "7c. Variable shear modulus"]]],
    ["Results", [
      ["official.html", "8. The official PTHA18 run"],
      ["data_provenance.html", "Data provenance"],
      ["report_guide.html", "Reading report.html"],
      ["validation.html", "Validation"]]],
    ["Your own runs", [
      ["zones.html", "Zones you can run"],
      ["recipes.html", "Recipes: one per zone"],
      ["customise.html", "Change a setting, re-run"],
      ["reference.html", "Options, files, columns"],
      ["troubleshooting.html", "Troubleshooting"]]],
    ["Code", [
      ["code_map.html", "Code map"],
      ["glossary.html", "Glossary"]]]
  ];

  var here = location.pathname.split("/").pop() || "index.html";
  var base = function (s) { return s.split("#")[0]; };
  var flat = [];
  GROUPS.forEach(function (g) { g[1].forEach(function (p) { flat.push(p); }); });

  // page list
  var side = document.getElementById("side");
  if (side) {
    var h = '<div class="brand"><a href="index.html">ptha_python</a></div>' +
            '<div class="sub">PTHA18 earthquake sources from public data</div>';
    GROUPS.forEach(function (g) {
      h += '<div class="grp">' + g[0] + '</div>';
      g[1].forEach(function (p) {
        h += '<a class="pg' + (base(p[0]) === here ? " active" : "") + '" href="' + p[0] + '">' + p[1] + '</a>';
      });
    });
    side.innerHTML = h;
  }

  // on this page
  var toc = document.getElementById("toc");
  var heads = Array.prototype.slice.call(document.querySelectorAll("main h2[id]"));
  if (toc && heads.length) {
    toc.innerHTML = '<div class="ttl">On this page</div>' + heads.map(function (x) {
      return '<a href="#' + x.id + '">' + x.textContent + '</a>';
    }).join("");
    var links = Array.prototype.slice.call(toc.querySelectorAll("a"));
    var mark = function () {
      var k = 0;
      heads.forEach(function (x, i) { if (x.getBoundingClientRect().top < 120) k = i; });
      links.forEach(function (a, i) { a.classList.toggle("on", i === k); });
    };
    document.addEventListener("scroll", mark, { passive: true });
    mark();
  }

  // previous / next
  var pager = document.getElementById("pager");
  var i = flat.map(function (p) { return base(p[0]); }).indexOf(here);
  if (pager && i >= 0) {
    var s = "";
    if (i > 0) s += '<a class="prev" href="' + flat[i - 1][0] + '"><span class="k">Previous</span>' + flat[i - 1][1] + '</a>';
    if (i < flat.length - 1) s += '<a class="next" href="' + flat[i + 1][0] + '"><span class="k">Next</span>' + flat[i + 1][1] + '</a>';
    pager.innerHTML = s;
  }

  // copy buttons
  document.querySelectorAll(".code").forEach(function (c) {
    if (c.querySelector(".copy")) return;
    var b = document.createElement("button");
    b.className = "copy"; b.textContent = "Copy";
    b.addEventListener("click", function () {
      var t = c.querySelector("pre").innerText;
      (navigator.clipboard ? navigator.clipboard.writeText(t) : Promise.reject()).then(function () {
        b.textContent = "Copied"; b.classList.add("done");
        setTimeout(function () { b.textContent = "Copy"; b.classList.remove("done"); }, 1400);
      }, function () {});
    });
    c.appendChild(b);
  });

  // every image opens at full size in a new tab
  document.querySelectorAll("main img").forEach(function (im) {
    if (im.closest("a")) return;
    var a = document.createElement("a");
    a.href = im.getAttribute("src"); a.target = "_blank"; a.title = "Open the full-size image";
    im.parentNode.insertBefore(a, im); a.appendChild(im);
    im.style.cursor = "zoom-in";
  });

  // tabs: <div class="tabs" data-group="x"><button data-pane="id">..</button></div> + <div class="pane" id="id">
  document.querySelectorAll(".tabs").forEach(function (t) {
    var btns = Array.prototype.slice.call(t.querySelectorAll("button[data-pane]"));
    btns.forEach(function (b) {
      b.addEventListener("click", function () {
        btns.forEach(function (o) {
          o.classList.toggle("on", o === b);
          var p = document.getElementById(o.getAttribute("data-pane"));
          if (p) p.classList.toggle("on", o === b);
        });
      });
    });
  });
})();
