module JS

using Gumbo, HTTP, Bonito
using Makie, WGLMakie
export FigureContainer

abstract type BackendType end
struct WGLMakieType <: BackendType end

function backend_type(::Module)
    if Makie.current_backend() === WGLMakie
        return WGLMakieType
    else
        return BackendType
    end
end

struct FigureContainer{T<:BackendType}
    fig::Figure
end
FigureContainer(fig::Figure) = FigureContainer{backend_type(Makie.current_backend())}(fig)

js_counter = 1
dir = @__DIR__

function separate_javascript!(element::HTMLElement)
    for (i, child) in enumerate(element.children)
        if (typeof(child) <: HTMLElement{:script}) &&
           haskey(child.attributes, "src") &&
           contains(child.attributes["src"], "data:application/javascript;base64,")
            src_raw = replace(
                child.attributes["src"],
                "data:application/javascript;base64," => "",
            )
            js_content = HTTP.base64decode(src_raw)
            global js_counter += 1
            js_filename = "./js/js-$(js_counter).js"
            js_write = "$dir/../_build/html/julia/js/js-$(js_counter).js"
            write(js_write, js_content)
            element.children[i] = HTMLElement{:script}(
                Vector{HTMLNode}[],
                element,
                Dict{String,String}("src" => js_filename, "type" => "module"),
            )
        end
        if typeof(child) <: HTMLElement
            separate_javascript!(child)
        end
    end
    return element
end

function separate_javascript(html_source::AbstractString)
    if !isdir("$dir/../_build/html/")
        mkdir("$dir/../_build/html/")
    end
    if !isdir("$dir/../_build/html/julia/")
        mkdir("$dir/../_build/html/julia/")
    end
    if !isdir("$dir/../_build/html/julia/js/")
        mkdir("$dir/../_build/html/julia/js/")
    end
    parsed_html = parsehtml(html_source)
    div_element = separate_javascript!(parsed_html.root.children[2].children[1])
end

function Base.display(fig::FigureContainer{<:BackendType})
    display(fig.fig)
end

function Base.display(fig::FigureContainer{WGLMakieType})
    html_str = string(DOM.div(fig.fig))
    html_source_new = separate_javascript(html_str)
    display("text/html", string(html_source_new))
end

end
