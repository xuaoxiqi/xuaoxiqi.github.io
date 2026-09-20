(() => {
    const page = document.querySelector('#articlePost.page-research');
    if (!page) return;

    const conferences = page.querySelector('[data-conf-section="conferences"]');
    const stats = page.querySelector('[data-conf-stats]');
    function updateConferenceStats() {
        if (!conferences || !stats) return;
        const featured = [...conferences.children].filter(item =>
            item.querySelector('.badge-invite, .badge-keynote')
        );
        stats.querySelector('[data-stat="featured"]').textContent = featured.length;
        stats.hidden = false;
    }
    updateConferenceStats();
    if (conferences) {
        new MutationObserver(updateConferenceStats).observe(conferences, {
            childList: true,
            subtree: true,
            attributes: true,
            attributeFilter: ['class']
        });
    }

    page.querySelectorAll('[data-start]').forEach(project => {
        if (new Date() < new Date(`${project.dataset.start}T00:00:00`)) return;
        const badge = project.querySelector('.status-badge');
        badge.textContent = 'Active';
        badge.classList.replace('badge-upcoming', 'badge-active');
    });

    const siteNav = document.querySelector('#navBarTop');
    const sectionNav = page.querySelector('.research-nav');
    if (!sectionNav) return;
    const links = [...sectionNav.querySelectorAll('a[href^="#"]')];
    const sections = links.map(link => document.querySelector(link.hash));
    let anchorOffset = 140;
    let framePending = false;

    function updateCurrentSection() {
        let current = 0;
        sections.forEach((section, index) => {
            if (section && section.getBoundingClientRect().top <= anchorOffset + 2) current = index;
        });
        links.forEach((link, index) => {
            if (index === current) link.setAttribute('aria-current', 'location');
            else link.removeAttribute('aria-current');
        });
        framePending = false;
    }

    // The site navigation wraps on small screens, so measure both sticky rows.
    function updateNavOffsets() {
        const siteHeight = siteNav ? siteNav.getBoundingClientRect().height : 0;
        anchorOffset = siteHeight + sectionNav.getBoundingClientRect().height + 20;
        page.style.setProperty('--site-nav-height', `${siteHeight}px`);
        page.style.setProperty('--research-anchor-offset', `${anchorOffset}px`);
        updateCurrentSection();
    }
    updateNavOffsets();
    if ('ResizeObserver' in window) {
        const observer = new ResizeObserver(updateNavOffsets);
        if (siteNav) observer.observe(siteNav);
        observer.observe(sectionNav);
    }
    window.addEventListener('resize', updateNavOffsets);
    window.addEventListener('scroll', () => {
        if (framePending) return;
        framePending = true;
        requestAnimationFrame(updateCurrentSection);
    }, { passive: true });
})();
