# The docs/ copy of a report adds site chrome after it is generated: the
# reports-home button and page navigation (57ff1b6, a4a431e) and the equation
# pop-up script. Compare the reports with that chrome taken out.
function strip_site_chrome(html::AbstractString)
    html = replace(html, r"<script defer src=\"[./]*assets/equation-glossary\.js\"></script>\n?" => "")
    html = replace(html, r"<nav [^>]*>.*?</nav>"s => "")
    html = replace(html, r"(?:<p>)?<a [^>]*href=\"(?:\.\./)+index\.html\"[^>]*>.*?</a>(?:</p>)?"s => "")
    replace(html, r"[ \t]*\n" => "\n")
end
