Comments 💬
-----------
.. usage in rst files:
    .. include:: /_static/comments_section.rst

.. usage in ipynb files: 
    1. you need to remove .. raw:: directive
    ```{include} /_static/comments_section.rst
    :heading-offset: 1
    ```
    1. with .. raw:: directive but you need to add a header.
    ```{eval-rst}
    .. include:: /_static/comments_section.rst
        :start-line: 3
    ```

.. comment:
    <iframe class="giscus-frame-light" title="Comments" scrolling="no" allow="clipboard-write">    </iframe>
    <div class="giscus-light">     </div>

Feel free to make a comment, ask a question, or share your opinion about this specific content. 
Please keep in mind the `Commenting Guidelines ⚖ </community/commenting_guidelines.html>`_.

.. raw:: html

    <script>
        let giscusTheme = localStorage.theme;
        let giscusAttributes = {
            "src": "https://giscus.app/client.js",
            "data-repo": "hiesabx/reservoirflow",
            "data-repo-id": "R_kgDOG3v41g",
            "data-category": "Documentation",
            "data-category-id": "DIC_kwDOG3v41s4CpxFL",
            "data-mapping": "pathname",
            "data-strict": "1",
            "data-reactions-enabled": "1",
            "data-emit-metadata": "0",
            "data-input-position": "bottom",
            "data-theme": giscusTheme,
            "data-lang": "en",
            "data-loading": "lazy",
            "crossorigin": "anonymous",
            "async": ""};
        let giscusScript = document.createElement("script");
        Object.entries(giscusAttributes).forEach(([key, value]) => giscusScript.setAttribute(key, value));
        document.getElementById("comments").appendChild(giscusScript);
    </script>

.. comments
    document.getElementById("comments").remove("giscus")
    document.getElementsByClassName("bd-content").appendChild
    document.getElementsByClassName("bd-article").appendChild

.. raw: html
    :class: only-light

    <script src="https://giscus.app/client.js"
            data-repo="hiesabx/reservoirflow_comments"
            data-repo-id="R_kgDOKTqNNg"
            data-category="General"
            data-category-id="DIC_kwDOKTqNNs4Cgs8l"
            data-mapping="pathname"
            data-strict="1"
            data-reactions-enabled="1"
            data-emit-metadata="0"
            data-input-position="bottom"
            data-theme="dark"
            data-lang="en"
            data-loading="lazy"
            crossorigin="anonymous"
            async
    >
    </script>

.. raw: html
    :class: only-dark

    <script src="https://giscus.app/client.js"
            data-repo="hiesabx/reservoirflow_comments"
            data-repo-id="R_kgDOKTqNNg"
            data-category="General"
            data-category-id="DIC_kwDOKTqNNs4Cgs8l"
            data-mapping="pathname"
            data-strict="1"
            data-reactions-enabled="1"
            data-emit-metadata="0"
            data-input-position="bottom"
            data-theme="light"
            data-lang="en"
            data-loading="lazy"
            crossorigin="anonymous"
            async
    >
    </script>

.. comment:
    .. raw:: html
        :class: only-dark

        <script 
            type="text/javascript"
            src="https://utteranc.es/client.js"
            async="async"
            repo="hiesabx/reservoirflow_comments"
            issue-term="pathname"
            theme="github-dark"
            label="comments 💬"
            crossorigin="anonymous"
        >
        </script>

    .. raw:: html
        :class: only-light

        <script 
            type="text/javascript"
            src="https://utteranc.es/client.js"
            async="async"
            repo="hiesabx/reservoirflow_comments"
            issue-term="pathname"
            theme="github-light"
            label="comments 💬"
            crossorigin="anonymous"
        >
        </script>
